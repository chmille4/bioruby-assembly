
require 'bio-assembly/read/ace'

module Bio
  class Assembly
    class Read 
      include Bio::Assembly::Read::Ace
      
      attr_accessor :seq, :name, :orientation, :from, :to, :clear_range_from, :clear_range_to
      def initialize(str="")
       @seq = Bio::Sequence::NA.new(str)
      end
      
      def ==(other_read)
         name == other_read.name
      end
      
      def num_bases
        seq.length
      end
      
      def from=(new_from)
        @from = new_from.to_i
      end
      
      def to=(new_to)
        @to = new_to.to_i
      end
      
      def clear_range_from=(new_clear_range_from)
        @clear_range_from = new_clear_range_from.to_i
      end
      
      def clear_range_to=(new_clear_range_to)
        @clear_range_to = new_clear_range_to.to_i
      end
      
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
      
      def <=>(other)
        unless other.kind_of?(Bio::Assembly::Read)
          raise "[Error] markers are not comparable"
        end
        if self.from == other.from
          # sort by to if froms are identical
          return self.to.<=>(other.to)
        else
          return self.from.<=>(other.from)
        end
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
      
    end
    
  end
end