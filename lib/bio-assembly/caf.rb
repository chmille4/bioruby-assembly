module Bio
  class Assembly
    class Caf < Bio::Assembly
      
      # register parser with superclass
      register_parser :caf

      def initialize(path)
         @file = File.new(path, 'r') 
      end
      # iterator that return one contig at a time
      def each_contig
        contig = Contig.new
        feature = Hash.new
        @file.each do |line|
          feature = parse_blocks(line,feature) # search the file for CAF blocks like DNA and Sequence
          if feature[:type] == :read and feature[:parsed]
            read = convert_to_read(feature)
            contig.add_read(read)
            feature = Hash.new
          elsif feature[:type] == :contig and feature[:parsed]
            contig = convert_to_contig(contig,feature)
            yield contig
            contig = Contig.new
            feature = Hash.new
          end
        end
      end

      class Contig < Bio::Assembly::Contig
      end 


      class Read < Bio::Assembly::Read
        attr_accessor :quality
      end
            
      private
      
      def parse_blocks(line,feat)
        keywords = line.split("\s")
        case keywords[0]
          when "DNA" then parse_dna(feat)
          when "Sequence" then parse_seq(feat,line)
        end  
        return feat
      end
      
      # parse DNA sequence and BaseQuality
      def parse_dna(feat)
        feat[:seq] = @file.gets("\n\n").tr("\n","")
        newline = @file.gets
        keywords = newline.split("\s")
        feat[:qual] = @file.gets("\n\n").tr("\n"," ").rstrip if keywords[0] == "BaseQuality"
        feat[:parsed] = true if feat[:type] == :contig
      end
      
      # parse Sequence information like Name, Clipping, Strand and Type
      def parse_seq(feat,line)
        feat[:name] = line.split(":")[1].tr("\s|\n","")
        sequence_block = @file.gets("\n\n")
        sequence_block.split("\n").each do |l|
          keywords = l.split("\s")
          case keywords[0]
            when "Clipping" then parse_clipping(feat,l)
            when "Strand" then parse_strand(feat,l)
            when "Assembled_from" then parse_af(feat,l)
            when "Is_read" then feat[:type] = :read
            when "Is_contig" then feat[:type] = :contig
          end
        end
        feat[:parsed] = true if feat[:type] == :read
      end
      
      # parse read coordinates for quality clipping
      def parse_clipping(feat,line)
        val = line.chomp.split("\s")
        feat[:clipping_start] = val[-2]
        feat[:clipping_end] = val[-1]
      end
      
      # parse sequence strand information
      def parse_strand(feat,line)
        feat[:orientation] = line.split("\s")[1].tr("\n","")
      end
      
      # parse Assembled_from lines in Contig. These lines also include read alignment positions within the contig
      def parse_af(feat,line)
        if feat[:af].nil?
          feat[:af] = [line]
        else
          feat[:af] << line
        end    
      end

      # convert a generic feature into a Caf::Read object
      def convert_to_read(feature)
        read = Read.new
        read.name = feature[:name]
        read.seq = feature[:seq]
        read.quality = feature[:qual]
        read.clear_range_from = feature[:clipping_start]
        read.clear_range_to = feature[:clipping_end]
        read.orientation = feature[:orientation]
        return read       
      end
      
      # convert a generic feature into a Caf::Contig object
      def convert_to_contig(contig,feature)
        contig.name = feature[:name]
        contig.seq = feature[:seq]
        contig.quality = feature[:qual]
        # assign reads ranges using Assembled_from lines in Contig
        feature[:af].each do |af|
          val = af.split("\s")
          contig.reads[val[-5]].from = val[-4]
          contig.reads[val[-5]].to = val[-3]
        end
        return contig       
      end

    end # end Caf
  end # end Assembly  
end # end Bio