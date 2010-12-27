
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
   
   # open contig class and write ace specific methods for contig objects         
   class Contig

      def to_ace
        ace = ""
        ace += ['CO', name, num_bases, num_reads, num_base_segments, orientation].join(' ') + "\n"
        ace += seq.to_s.gsub(Regexp.new(".{1,50}"), "\\0\n") + "\n"
        ace += "BQ\n"
        last_stop = quality.size - 1
        (quality.size/50+1).times do |i|
          start = i * 50
          stop = (i+1) * 50 - 1
          stop = last_stop if stop > last_stop
          ace += ' ' + quality[start..stop].join(' ')  + "\n"
        end
        ace += "\n"

        # holds BS data for reads
        bs_str = ""
        # holds RD, QA, and DS data for reads 
        rest_str = ""
        @reads.values.sort.each do |read|
          ace += read.to_ace_af
          bs_str += read.to_ace_bs
          rest_str += read.to_ace_rest
        end

        # compile data in correct order
        ace += bs_str
        ace += "\n"
        ace += rest_str
        ace
      end
         
   end # => end Contig class
   
   # open Read class to add ace specific methods for read objects
   class Read

     attr_accessor :base_sequences
  
     def to_ace
       ace += ""
       # holds BS data for reads
       bs_str = ""
       # holds RD, QA, and DS data for reads 
       rest_str = ""
       ace += to_ace_af
       bs_str += to_ace_bs
       rest_str = to_ace_rest

       # compile data in correct order
       ace += bs_str
       ace += "\n"
       ace += rest_str
       ace
     end

     def to_ace_bs
       bs_str = ""
       unless base_sequences.nil?
         base_sequences.each do |bs|
           bs_str += ['BS', bs.from, bs.to, bs.read_name].join(' ') + "\n"
         end
       end
       bs_str
     end

     def to_ace_af
       ['AF', name, orientation, from].join(' ') + "\n"
     end

     def to_ace_rest
       rest_str = ""
       rest_str += ['RD', name, num_bases, 0, 0].join(' ') + "\n"
       rest_str += seq.to_s.gsub(Regexp.new(".{1,50}"), "\\0\n")  + "\n"
       rest_str += ['QA', clear_range_from, clear_range_to, clear_range_from, clear_range_to].join(' ') + "\n"
       rest_str += ['DS', 'CHROMAT_FILE:', name, 'PHD_FILE:', "#{name}.phd.1", 'TIME:', Time.now].join(' ') + "\n"
       rest_str
     end
  
     def add_base_sequence(from, to, read_name)
       @base_sequences = Array.new if @base_sequences.nil?
       @base_sequences.push BaseSequence.new(from, to, read_name)
     end
  
     class BaseSequence
        attr_accessor :from, :to, :read_name
    
        def initialize(from, to, read_name)
           @from = from
           @to = to
           @read_name = read_name
        end
    
        def <=>(other)
           unless other.kind_of?(Bio::Assembly::Read::BaseSequence)
              raise "[Error] markers are not comparable"
           end
           if self.from == other.from
              # sort by to if froms are identical
              return self.to.<=>(other.to)
           else
              return self.from.<=>(other.from)
           end
        end    
  
      end # => end BaseSequence Class

   end # => end Read Class
   

end # => end class Assembly
end # => end module Bio