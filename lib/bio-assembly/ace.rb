
module Bio
class Assembly

class Ace < Bio::Assembly

   # register parser with superclass
   register_parser :ace

   def initialize(path)
      @file = File.new(path, 'r') 
      @contigs = Array.new
      parse_as
    end

   def each_contig
     # check if file is already parsed
     if @total_num_contigs.to_i == @contigs.size
       @contigs.each{ |contig| yield contig }
     else
       each_identifier do |identifier, attrs|
         next unless identifier == 'CO'
         contig = parse_contig(attrs)
         @contigs.push contig
         yield(contig)   
       end
     end
   end

   def to_ace
     ace = ""
     ace += "AS " + num_contigs.to_s + " " + num_reads.to_s + "\n\n"
     each_contig { |contig| ace += contig.to_ace + "\n" }
     ace
   end

   private
   def parse_contig(attrs)
     contig = Bio::Assembly::Contig.new
     contig.name, base_num, @num_reads, base_segments_num, contig.orientation = attrs.split(" ")
     # keep track of the number of RD identifiers parsed
     @num_rds_parsed = 0

     # get sequence
     seq = @file.gets("\n\n").tr(" \r\n", "")
     contig.seq = seq

     # loop through identifiers (e.g AF, RD, etc)
     each_identifier do |identifier, attrs|    
       case identifier
         when "BQ" then parse_bq(contig)
         when "AF" then parse_af(contig, attrs)
         when "BS" then parse_bs(contig, attrs)
         when "RD" then parse_rd(contig, attrs); break if @num_rds_parsed == @num_reads.to_i
         when "WR" then parse_wr(contig, attrs)
         when "RT" then parse_rt(contig, attrs)
         when "CT" then parse_ct(contig, attrs)
         when "WA" then parse_wa(contig, attrs)
       end
     end

    contig
   end

   # Finds the next_identifier
   def each_identifier
     @file.each do |line|
       next if line !~ /^[ABCDQRW][ADFOQRST][\s\n].*/
       yield(line[0..1], line[3..-1])
     end
   end

   # parse assembly meta data
   def parse_as
     line = @file.gets
     identifier, @total_num_contigs, total_num_reads = line.split(" ")
   end

   # parse contig sequence quality data
   def parse_bq(contig)
     contig.quality = @file.gets("\n\n").tr("\r\n", "").gsub(/^\s/, "").split(' ')
   end

   # parse read meta data
   def parse_af(contig, attrs)
     read = Bio::Assembly::Read.new
     read.name , read.orientation, read.from = attrs.split(" ")
     contig.add_read read
   end

   # parse base sequence data
   def parse_bs(contig, attrs)
     from, to, read_name = attrs.split(" ")
     read = contig.find_read_by_name( read_name )
     read.add_base_sequence(from, to, read_name)
   end

   # parse read sequence and position data
   def parse_rd(contig, attrs)
     # increment counter
     @num_rds_parsed += 1

     # parse read
     read_name, num_padded_bases, num_read_infos, num_read_tags = attrs.split(" ") 
     seq = @file.gets("\n\n").tr( " \r\n", "")

     # get read with matching name
     read = contig.find_read_by_name( read_name )
     read.seq = seq
     read.to = read.from.to_i + read.seq.length
     # set read.to to contig length if read runs off contig
     read.to = contig.seq.length if read.to > contig.seq.length

     # if present parse QA and DS associated with this read
     each_identifier do |identifier, attrs|
       case identifier
         when "QA" then parse_qa(read, attrs)
         when "DS" then parse_ds(read, attrs); break
       end
     end

   end

   # parse a read's clear ranges (the part of the read that contributes to the contig)
   def parse_qa(read, attrs)
     start, stop, clear_range_from, clear_range_to = attrs.split(" ")
     read.clear_range_from = clear_range_from
     read.clear_range_to = clear_range_to
   end

   # parse file data - ignored
   def parse_ds(read, attrs)
   end

   # parse run meta data - ignored
   def parse_wa(contig, attrs)
   end

   # parse run meta data - ignored
   def parse_ct(contig, attrs)
   end

end # => end class Ace

end # => end module Assembly
end # => end module Bio