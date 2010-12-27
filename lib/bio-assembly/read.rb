
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
      
    end
    
  end
end