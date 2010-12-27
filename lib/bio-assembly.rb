require 'bio/sequence' 
require 'bio-assembly/contig'
require 'bio-assembly/read'

module Bio

   class Assembly
      attr_accessor :contigs
      
      @@formats = { }
    
      def self.create(path, format)
         streamer = @@formats[format]
         if streamer
            streamer.new(path)
         else
            raise "Format type '#{format}' is not supported"
         end
      end
      
      def self.register_parser name
         @@formats[name] = self
      end
      
      def contigs
        # use each_contig to stream large files
        parse_whole_file if @contigs.empty?
        @contigs
      end
      
      def each_contig
         # implemented by each format subclass
      end
      
      private
      
      def num_contigs
        contigs.size
      end

      def num_reads
        read_num = 0
        each_contig { |contig| read_num += contig.num_reads }
        read_num
      end

      def parse_whole_file
        each_contig { |x| 1 }
      end
    
   end
   
end

require 'bio-assembly/ace'