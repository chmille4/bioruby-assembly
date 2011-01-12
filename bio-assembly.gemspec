# Generated by jeweler
# DO NOT EDIT THIS FILE DIRECTLY
# Instead, edit Jeweler::Tasks in Rakefile, and run 'rake gemspec'
# -*- encoding: utf-8 -*-

Gem::Specification.new do |s|
  s.name = %q{bio-assembly}
  s.version = "0.1.0"

  s.required_rubygems_version = Gem::Requirement.new(">= 0") if s.respond_to? :required_rubygems_version=
  s.authors = ["Chase Miller", "Francesco Strozzi"]
  s.date = %q{2011-01-12}
  s.description = %q{bioruby plugin to parse, write, and manipulate assembly data}
  s.email = %q{chmille4@gmail.com}
  s.extra_rdoc_files = [
    "LICENSE.txt",
    "README.rdoc"
  ]
  s.files = [
    ".document",
    "Gemfile",
    "Gemfile.lock",
    "LICENSE.txt",
    "README.rdoc",
    "Rakefile",
    "VERSION",
    "bio-assembly.gemspec",
    "data/example.caf",
    "data/example1.ace",
    "lib/bio-assembly.rb",
    "lib/bio-assembly/ace.rb",
    "lib/bio-assembly/caf.rb",
    "lib/bio-assembly/contig.rb",
    "lib/bio-assembly/maf.rb",
    "lib/bio-assembly/read.rb",
    "test/helper.rb",
    "test/test_bio-assembly-ace.rb",
    "test/test_bio-assembly-caf.rb"
  ]
  s.homepage = %q{http://github.com/chmille4/bioruby-assembly}
  s.licenses = ["MIT"]
  s.require_paths = ["lib"]
  s.rubygems_version = %q{1.3.7}
  s.summary = %q{BioRuby Assembly plugin}
  s.test_files = [
    "test/helper.rb",
    "test/test_bio-assembly-ace.rb",
    "test/test_bio-assembly-caf.rb"
  ]

  if s.respond_to? :specification_version then
    current_version = Gem::Specification::CURRENT_SPECIFICATION_VERSION
    s.specification_version = 3

    if Gem::Version.new(Gem::VERSION) >= Gem::Version.new('1.2.0') then
      s.add_runtime_dependency(%q<bio>, [">= 1.4.1"])
      s.add_development_dependency(%q<shoulda>, [">= 0"])
      s.add_development_dependency(%q<bundler>, ["~> 1.0.0"])
      s.add_development_dependency(%q<jeweler>, ["~> 1.5.2"])
      s.add_development_dependency(%q<rcov>, [">= 0"])
      s.add_development_dependency(%q<bio>, [">= 1.4.1"])
    else
      s.add_dependency(%q<bio>, [">= 1.4.1"])
      s.add_dependency(%q<shoulda>, [">= 0"])
      s.add_dependency(%q<bundler>, ["~> 1.0.0"])
      s.add_dependency(%q<jeweler>, ["~> 1.5.2"])
      s.add_dependency(%q<rcov>, [">= 0"])
      s.add_dependency(%q<bio>, [">= 1.4.1"])
    end
  else
    s.add_dependency(%q<bio>, [">= 1.4.1"])
    s.add_dependency(%q<shoulda>, [">= 0"])
    s.add_dependency(%q<bundler>, ["~> 1.0.0"])
    s.add_dependency(%q<jeweler>, ["~> 1.5.2"])
    s.add_dependency(%q<rcov>, [">= 0"])
    s.add_dependency(%q<bio>, [">= 1.4.1"])
  end
end

