require 'helper'

class TestBioAssemblyCaf < Test::Unit::TestCase
  
  def setup
    path = File.join('data','example.caf')
    @caf = Bio::Assembly.open(path,:caf)
    @contigs = []
    @caf.each_contig {|c| @contigs << c}
  end
  
  def test_contigs
    assert_equal(2,@contigs.size)
    assert_equal('Contig1',@contigs[0].name)
    assert_equal('Contig2',@contigs[1].name)
    assert_instance_of(Bio::Assembly::Caf::Contig,@contigs[0])
  end
  
  def test_read_per_contig
    assert_equal(21,@contigs[0].reads.size)
    assert_equal(21,@contigs[0].num_reads)
  end
  
  def test_contig_seq
    seq = "TTCAAGCGATTCTCCTGCCTCAGCCTCCCGAGTAGCTGGGATTATANACTGTGCGTGCGCCACCATGCCTGGCTAATTTTTGTATTTTTAGTAGGGATGGGGTTTCACCATGTTGGCCAAGCTGGTCTCGAGCTCCTGACCTAGGATTACAGGCCTAAGCCACCGCACCCGGCATGATGGGTCTTTATTCTTCAAAGCAGGAGGAAGGGATCCTAGAAAAACAGAGACAAGGCCAAACATGGTAGCTCACACCTGTAATNNNNNCACTTTGGGAGGCCAGTGCGGGTGAATCACGANGTCAGGAGTTCAAGACCAGCCTGGCCAACATGGTGAANCACCGTCTCTACTAAAATAAAAAGAAATTAGCTGGGTGTCGTGGCAGGTGCNTGTAATCCCAGCCACTTGGGAAGCTGAGGCAGGAGAATCGCTTGAACCCGGGAGGTGGAGGTTGCAGTGAGCCGAGATCACGCGACTGCACTCCAGCCCAACCAATAGTGTGAGACTCTGTCTCGAAAAAAAAAAAGCAGAGACAAGACNACTAGTACAGTACTTACAGGGTTATTATGATGATTAAATGAGAGAATAGCTGTGAGGTGATTGATATAGTGCTGTGCTTAATACAAACTATCATTTTATTATACGGGTTGAGTGTNTCTAATCTGAAAATCCAAAATTAGAAATGCTCTACAGTCTGAAACTTTTTTGAGCACCGACCTAATGTTCAAAGGAAGTGCTTATTGGAGCATTATGGGTTGTTAGATTTTTGGGTTGGGAATATTCAACCAGTAAGTACTATAAAATGCAAATATTCCAAAAAAAATCTGAAATCTGAAACATTTCTGGTCCTAAGCAAGCATTTTGCAAAGGGATACGCAACCTGTAGTACGTTCTTTATCATTGTTTTAAGTAGTTAATATATTGTGGTACAGATTCTGAGGTGGTATAGCAAATTCGATTGTATTATTAAAAAGCATATTTATATTTTGAGAGCTTGCTTAGGATTATTGGAGAGAATAAAACAGTGAAGCTTTGGTGTTATGAGGGAATTTTAGATAGAAAAGTGCAGTTTTTCAGTTCATGCTCTTTCATTTTTTACTCCCTCAGGTTAAAGCTNGAAGCTCAACAAAGATATAGTGATCTCTGTGGGCATTTATAATCTGGTCCAGAAGGCTCTNNANNCNNNTCCNNNNNNNCTNNANNNNNNNACAAATGAACCAGTGAAAACCAAGACCCGGACCTTTAATACAAGTACAGGCGGTTTGCTTCTGCCTAGCGATACCAAGAGGTCTCAGGTAGGTAGAGATGCCTTTTGTTGTTGTTGTTTTTGAGACAGGGTCTCATTGTGTCGCCCAGGCTGGAGTGCAGTGGGGCGAACATGACTCGCTACAGCCTTGACCTCCTGGACTCAAGCGATCCTTCTGTCTCAGCCTCCCAAGTAGCTGGGATCACAGGCATGTGACATCACACCCAGCTAATTTATTTATTTATTTATTTTTTAAGAGACTGGATCGACTGGGCACAGTGGCTCATGCCTGTAATCCCANCACTTTGGGAGGCCGAGGCAGGTGGATTACCGAGGTCAGGAGTTCAAGACCAGCCTGACCAACATGGAGAAACCCCATCTCTACTAAAAATACAAAATGAGCTGGGCATGGTGGTGCATGCCTGTAATCC"
    assert_equal(seq,@contigs[0].seq.to_s.upcase)
    assert_instance_of(Bio::Sequence::NA,@contigs[0].seq)
  end
  
  def test_contig_qual
    qual = "4 4 8 4 6 10 13 21 24 25 33 33 33 30 27 21 15 27 19 30 30 33 33 30 21 21 10 9 17 11 27 27 37 38 33 35 35 35 35 44 45 45 38 37 37 45 45 45 43 43 43 45 45 45 45 45 45 45 21 21 23 30 30 34 37 37 38 45 32 45 37 41 37 45 30 45 45 34 34 34 32 29 22 32 32 45 45 45 33 28 28 37 37 34 34 34 35 34 34 34 34 34 34 37 37 40 41 37 34 37 37 37 32 24 22 27 29 25 27 20 21 21 21 27 45 41 40 40 40 42 34 34 37 41 45 51 45 45 37 30 37 41 41 37 36 28 30 30 30 22 33 33 35 33 33 41 51 51 45 39 39 39 30 28 33 34 34 39 34 41 34 34 41 33 33 39 39 39 33 33 33 30 33 33 33 30 29 19 19 24 25 32 33 45 36 36 36 36 41 41 30 30 30 37 37 43 51 51 51 51 45 37 37 30 37 37 37 51 51 51 51 51 37 37 28 28 10 10 10 13 10 10 10 9 9 9 21 21 28 37 37 37 33 33 33 33 33 33 34 34 34 33 33 33 33 21 25 25 22 22 29 29 33 33 33 33 31 31 27 28 17 23 23 28 26 24 32 10 10 10 15 15 26 32 32 37 45 45 32 45 33 26 24 32 37 37 35 34 37 35 37 34 37 37 40 37 37 45 45 38 38 45 45 40 37 37 37 34 30 30 28 28 28 30 37 37 34 34 34 34 34 33 31 31 25 21 19 22 21 25 30 32 32 30 22 22 25 28 24 25 26 17 17 23 27 27 33 31 31 34 38 38 38 38 34 26 17 15 8 8 11 16 25 31 32 33 29 26 32 29 33 34 34 34 36 33 34 34 31 27 31 34 38 45 45 36 37 37 37 36 36 38 34 32 28 25 22 22 21 17 18 29 29 23 23 22 27 25 20 18 11 11 15 15 16 23 21 17 19 24 24 24 31 31 31 33 33 33 31 33 33 29 29 29 29 29 31 31 31 23 24 19 12 10 16 10 10 20 10 10 13 15 30 23 29 23 28 18 10 10 16 21 18 19 19 24 24 12 11 9 9 10 21 23 31 31 33 28 28 14 17 17 28 21 30 24 30 28 22 26 23 19 10 10 12 23 22 23 19 12 10 9 9 10 19 24 29 30 34 34 34 34 34 29 25 25 32 26 31 31 32 20 17 15 12 12 4 4 4 17 23 33 30 24 18 13 18 16 12 16 8 9 9 19 15 15 4 4 4 7 11 6 6 6 7 8 9 17 12 13 19 23 25 13 13 13 19 21 29 32 30 26 20 20 12 9 9 8 9 8 8 15 24 17 15 8 8 9 17 16 4 4 4 8 8 13 15 23 14 9 8 8 10 10 17 23 21 17 12 13 22 20 19 15 11 11 9 8 8 8 8 8 8 10 9 9 7 8 8 10 12 10 10 12 4 4 4 4 4 4 10 10 10 10 10 9 8 8 7 8 9 11 9 9 4"
    assert_equal(qual,@contigs[1].quality)
  end
  
  def test_read
    read = nil
    assert_nothing_raised do 
      read = @contigs[0].find_read_by_name("22ak93c2.r1t")
    end
    seq = "GTCGCNCATAAGATTACGAGATCTCGAGCTCGGTACCCTTCAAGCGATTCTCCTGCCTCAGCCTCCCGAGTAGCTGGGATTATAGACTGTGCGTGCGCCACCATGCCTGGCTAATTTTTGTATTTTTAGTAGGGATGGGGTTTCACCATGTTGGCCAAGCTGGTCTCGAGCTCCTGACCTAGGATTACAGGCCTAAGCCACCGCACCCGGCATGATGGGTCTTTATTCTTCAAAGCAGGAGGAAGGGATCCTAGAAAAACAGAGACAAGGCCAAACATGGTAGCTCACACCTGTAATCCCANCACTTTGGGAGGCCAGTGCGGGTGAATCACGAAGTCAGGAGTTCAAGACCACCCTGGCCAACATGGTGTAACACCTGCTCTCCTAAAATTAAACAAAATTTCATGGTTTGCGTGGGCCGTCTTGTCTCATCACTTCACTCCTGAGGGCCGGCGCCGGAAAGATATCTTGATCTGCGGCGCTCCGACCGTTTTCTTTAAACCTTACAACTCCCGACCTCCTCGCCTATCCTCCCTAAATCCTCGCCAGGCTCGCCTGCTTCAGCCACTCTTTCCTTCGCACCCTCCCCTCTCTTCAATATACTTCACCCGCCCATCCTTCACGCCGGCACGTATCCAATCTCTTCTTATCTTTCCGTATCCAANTCCCTTCTCCCTCTGCCGCGACCTTCGCCATCCCTCTGCGCGTCCTCTTCC"
    assert_equal(seq,read.seq.to_s.upcase)
    assert_instance_of(Bio::Sequence::NA,read.seq)
    qual = "4 4 8 4 4 4 4 4 4 4 4 4 6 8 17 21 14 7 6 6 6 7 7 6 8 14 16 21 15 20 20 24 26 21 18 18 14 14 19 23 10 8 8 15 20 16 29 26 34 29 39 29 31 29 31 34 32 27 27 25 19 19 24 31 33 36 34 34 34 26 27 22 32 32 36 28 28 15 15 15 28 28 34 30 12 12 22 27 31 31 31 31 31 23 24 27 21 24 24 29 27 27 27 34 34 36 38 38 38 36 36 40 36 37 38 45 45 36 34 33 31 31 34 34 33 33 28 28 27 23 24 11 11 10 10 18 25 21 20 17 17 17 20 15 24 18 24 26 23 23 18 20 25 23 30 30 30 33 33 37 37 32 37 37 32 45 35 37 37 37 40 36 49 49 36 36 34 33 20 15 9 9 8 7 12 22 21 28 28 30 33 36 36 36 34 31 31 25 31 28 26 26 24 20 17 9 11 8 9 10 23 23 31 23 23 15 9 9 15 33 26 33 33 31 25 25 22 31 24 23 12 10 12 11 9 8 9 7 7 8 8 9 18 12 9 9 18 20 26 31 21 21 9 8 8 11 13 21 21 23 15 15 15 15 15 17 17 9 7 9 19 20 21 21 25 25 25 25 25 23 23 9 9 9 21 16 24 24 24 24 26 33 33 33 31 31 27 15 17 7 4 4 4 16 20 27 33 34 34 23 15 14 8 9 6 6 9 1 14 16 8 11 15 23 25 34 36 31 33 16 16 6 6 6 9 8 14 9 11 11 9 13 13 10 8 10 9 9 7 8 20 20 20 14 14 10 10 10 16 8 8 6 6 8 9 10 7 8 8 8 8 6 8 6 6 8 12 9 8 7 13 10 8 8 9 8 8 8 9 6 6 6 6 8 6 6 6 6 6 10 11 10 12 12 10 7 6 6 7 6 6 7 8 7 6 6 7 6 8 6 6 8 8 8 8 6 6 6 8 6 6 8 6 6 11 15 9 9 9 9 9 9 9 9 9 9 8 9 9 6 6 6 7 9 6 6 6 7 7 7 8 9 7 7 9 9 11 8 11 11 11 9 10 9 8 8 8 10 6 6 9 8 6 6 6 6 8 8 9 11 9 8 8 8 8 8 8 8 9 9 9 10 10 8 6 6 8 6 6 10 11 10 6 6 8 6 6 6 6 8 11 8 9 8 8 8 8 8 9 6 6 9 8 8 8 6 6 13 15 12 12 8 8 6 6 8 9 8 8 9 8 8 6 6 6 8 8 6 6 6 8 8 8 6 6 8 6 6 8 6 6 9 8 8 6 6 8 8 6 6 8 6 7 8 9 10 11 11 10 9 9 8 10 9 10 10 8 9 13 8 8 8 13 13 11 7 7 7 11 8 8 10 11 10 9 9 14 7 7 10 10 10 8 8 8 8 6 7 8 8 6 6 6 6 6 6 6 6 6 8 9 10 8 8 6 6 8 8 8 8 8 8 8 8 8 8 6 6 8 8 9 10 10 4 4 4 6 7 6 6 6 6 10 8 8 9 15 11 6 6 8 8 8 8 8 8 9 9 8 9 7 7 8 6 6 8 13 10 8 8 7 6 8 6 6 8 6 7 8 8 8 8 8 8 8 4"
    assert_equal(qual,read.quality)
    assert_equal(39,read.clear_range_from)
    assert_equal(331,read.clear_range_to)
    assert_equal(1,read.from)
    assert_equal(293,read.to)
    assert_equal("Reverse",read.orientation)
    assert_instance_of(Bio::Assembly::Caf::Read,read)
  end
  
  
end