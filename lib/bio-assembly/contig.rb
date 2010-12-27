module Bio
  class Assembly
    
    class Contig 
      attr_accessor :seq, :orientation, :quality, :to, :from, :name, :reads
      alias consensus_seq seq

      def initialize(str="")
        @reads = Hash.new
        @seq = Bio::Sequence::NA.new(str)
        # counter for RD identifier
        @rds_parsed = 0
      end

      def find_read_by_name(name)
        @reads[name]
      end
      
      def find_reads_in_range(clear_range_from, clear_range_to)
        reads_in_range = Array.new
        each_read do |read|
          
          # Read starts in region
          if read.from+read.clear_range_from > clear_range_from and read.from+read.clear_range_from < clear_range_to
             reads_in_range.push read
          # Read ends in region
          elsif read.to+read.clear_range_to < clear_range_to and read.to+read.clear_range_to > clear_range_from
             reads_in_range.push read
          # Read encompasses region
          elsif read.from+read.clear_range_from < clear_range_from and read.to+read.clear_range_to > clear_range_to
             reads_in_range.push read
          end
          
        end
        reads_in_range;
      end
      
      def add_read(read)
        # TODO do some checks for pos location
        @reads[read.name] = read
      end
  
      def each_read
        @reads.each_value { |read| yield read }
      end
      
      def num_reads
        @reads.size
      end
      
      def num_bases
        seq.length
      end
      
      def num_base_segments
        num_base_sequences = 0
        each_read do |read|
          num_base_sequences += read.base_sequences.size unless read.base_sequences.nil?
        end
        num_base_sequences
      end
      
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

    end
    
  end
end