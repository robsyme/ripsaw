#!/usr/bin/env ruby

require 'bio'
require 'pp'

def number_with_delimiter(number, delimiter=",", separator=".")
  begin
    parts = number.to_s.split('.')
    parts[0].gsub!(/(\d)(?=(\d\d\d)+(?!\d))/, "\\1#{delimiter}")
    parts.join separator
  rescue
    number
  end
end

filename = ARGV.shift
factory = Bio::Blast.local('blastn', filename, '-e 1e-10')

db = Hash[Bio::FlatFile.open(filename).map{|entry| [entry.entry_id, entry.naseq]}]

genome_size = 0
db.each do |entry_id, seq|
  seq.composition.each do |nuc, count|
    genome_size += count unless nuc == "n"
  end
end

$stderr.puts "Genome size: #{number_with_delimiter(genome_size)} bp"
`formatdb -p F -i #{filename}` unless ['.nhr', '.nin', '.nsq'].all?{|ext| File.exist?(filename + ext)}

i = 0
cursor = 0.0
ARGF.each do |line|
  seqid, start, stop, _ = line.split("\t")
  start = start.to_i
  stop = stop.to_i
  cursor += (stop - start + 1)
  $stderr.puts "%s:%d-%d (%s bp - %.2f%%)" % [seqid, start, stop, number_with_delimiter(stop - start + 1), cursor / genome_size * 100]
  fragment = db[seqid].subseq(start+1, stop.to_i)
  report = factory.query(fragment)
  seqs = [fragment]
  report.iterations.first.hits.each do |hit|
    hit.hsps.each do |hsp|
      if hsp.hit_frame > 0
        seqs << Bio::Sequence::NA.new("-" * hsp.query_from) + hsp.hseq
      else
        seqs << Bio::Sequence::NA.new("-" * hsp.query_to) + Bio::Sequence::NA.new(hsp.hseq).complement
      end
    end
  end

  if seqs.length > 1
    aln = Bio::Alignment.new(seqs).normalize
    sites = Enumerator.new{|yielder| aln.each_site{|site| yielder << site if site[0] != "-"}}
    sites.each_with_index do |site, j|
      depth = site.count{|nuc| nuc != "-" && nuc != "n"}
      next unless depth > 10
      counts = {"a" => 0, "c" => 0, "g" => 0, "t" => 0, "n" => 0, "-" => 0}
      site.each{|e| counts[e] += 1}
      ct = (counts["c"] + counts["t"]) / depth.to_f
      ag = (counts["a"] + counts["g"]) / depth.to_f
      ac = (counts["a"] + counts["c"]) / depth.to_f
      gt = (counts["g"] + counts["t"]) / depth.to_f
      bias = [ct,ag].max / [ac,gt].max
      if bias > 1
        if ct > ag
          # Change 't' to 'c'
          if db[seqid][start+j] == 't'
            db[seqid][start+j] = 'c'
            puts [seqid, start+j, start+j+1, "t->c", "%.4f" % bias].join("\t")
          end
        else
          # Change 'g' to 'a'
          if db[seqid][start+j] == 'a'
            db[seqid][start+j] = 'g'
            puts [seqid, start+j, start+j+1, "g->a", "%.4f" % bias].join("\t")
          end
        end
      end
    end
  end
end

File.open(filename + ".deripped", 'w') do |file|
  db.each do |entry_id, seq|
    file.puts seq.to_fasta(entry_id, 80)
  end
end
