CREATE TABLE "gly1" ("sequence_id" TEXT,
  "sequence" TEXT,
  FOREIGN KEY (sequence_id)
       REFERENCES summary (sequence_accession) 
       ON UPDATE CASCADE 
       ON DELETE CASCADE
);
CREATE INDEX "ix_gly1_sequence_id" ON "gly1" ("sequence_id");

CREATE TABLE "gly2" ("sequence_id" TEXT,
  "sequence" TEXT,
  FOREIGN KEY (sequence_id)
       REFERENCES summary (sequence_accession) 
       ON UPDATE CASCADE 
       ON DELETE CASCADE       
);
CREATE INDEX "ix_gly2_sequence_id" ON "gly2" ("sequence_id");

CREATE TABLE "gly3" ("sequence_id" TEXT,
  "sequence" TEXT,
  FOREIGN KEY (sequence_id)
       REFERENCES summary (sequence_accession) 
       ON UPDATE CASCADE 
       ON DELETE CASCADE       
);
CREATE INDEX "ix_gly3_sequence_id" ON "gly3" ("sequence_id");

CREATE TABLE "glyx3" ("sequence_id" TEXT,
  "sequence" TEXT,
  FOREIGN KEY (sequence_id)
       REFERENCES summary (sequence_accession) 
       ON UPDATE CASCADE 
       ON DELETE CASCADE       
);
CREATE INDEX "ix_glyx3_sequence_id" ON "glyx3" ("sequence_id");

CREATE TABLE "nterdb" ("nter" TEXT,
  "sequence_id" TEXT,
  "sequence" TEXT,
  FOREIGN KEY (sequence_id)
       REFERENCES summary (sequence_accession) 
       ON UPDATE CASCADE 
       ON DELETE CASCADE       
);
CREATE INDEX "ix_nterdb_sequence_id" ON "nterdb" ("sequence_id");

CREATE TABLE "summary" ("sequence_accession" TEXT,
  "sequence_src" TEXT,
  "flag" TEXT,
  "nter" TEXT,
  "nter_neighbor" TEXT,
  "cter" TEXT,
  "sequence" TEXT
);
CREATE INDEX "ix_summary_sequence_id" ON "summary" ("sequence_accession");

CREATE TABLE "features" ("sequence_id" TEXT,
  "sequence_src" TEXT,
  "feature_type" TEXT,
  "feature_start" INTEGER,
  "feature_end" INTEGER,
  "feature_id" TEXT,
  "pident" REAL,
  "coverage" REAL,
  "e-value" REAL,
  "feature_src" TEXT,
  "feature_target_len" INTEGER,
  "feature_seq" TEXT,
  FOREIGN KEY (sequence_id)
       REFERENCES summary (sequence_accession) 
       ON UPDATE CASCADE 
       ON DELETE CASCADE       
);
CREATE INDEX "ix_features_sequence_id"ON "features" ("sequence_id");

CREATE TABLE "hits" ("sequence_id" TEXT,
  "sequence_src" TEXT,
  "hit_target_len" INTEGER,
  "hit_start" INTEGER,
  "hit_stop" INTEGER,
  "hit_pident" REAL,
  "hit_e_value" REAL,
  "hit_coverage" REAL,
  "hit_method" TEXT,
  "hit_src" TEXT,
  FOREIGN KEY (sequence_id)
       REFERENCES summary (sequence_accession) 
       ON UPDATE CASCADE 
       ON DELETE CASCADE       
);
CREATE INDEX "ix_hits_sequence_id" ON "hits" ("sequence_id");

CREATE TABLE "ccya" ("sequence_id" TEXT,
  "ccyA_genomic_region" TEXT,
  "ccyA_start" INTEGER,
  "ccyA_stop" INTEGER,
  "ccyA_frame" TEXT,
  "ccyA_partial" TEXT,
  "ccyA_pseudo" TEXT,
  "ccyA_src" TEXT,
  "ccyA_seq" TEXT,
  FOREIGN KEY (sequence_id)
       REFERENCES summary (sequence_accession) 
       ON UPDATE CASCADE 
       ON DELETE CASCADE       
);
CREATE INDEX "ix_ccya_sequence_id"ON "ccya" ("sequence_id");

CREATE TABLE "harley" (
  "Accession" TEXT,
  "Date" TEXT
);
CREATE INDEX "ix_harley_accession" ON "harley" ("Accession");

CREATE TABLE "genomes" (
"Accession" TEXT,
  "Assembly name" TEXT,
  "Submitter" TEXT,
  "Submission date" TEXT,
  "Isolate" REAL,
  "TaxID" INTEGER,
  "Organism" TEXT,
  "Biosample" TEXT,
  "Isolation source" TEXT,
  "Environment (biome)" TEXT,
  "Geographic location" TEXT,
  "Culture collection" TEXT,
  "Collection date" TEXT,
  "Sample type" TEXT
);
CREATE INDEX "ix_genomes_accession" ON "genomes" ("Accession");

CREATE TABLE "gtdbtk" (
"user_genome" TEXT,
  "classification" TEXT,
  "fastani_reference" TEXT,
  "fastani_reference_radius" INTEGER,
  "fastani_taxonomy" TEXT,
  "fastani_ani" REAL,
  "fastani_af" REAL,
  "closest_placement_reference" TEXT,
  "closest_placement_radius" REAL,
  "closest_placement_taxonomy" TEXT,
  "closest_placement_ani" REAL,
  "closest_placement_af" REAL,
  "pplacer_taxonomy" TEXT,
  "classification_method" TEXT,
  "note" TEXT,
  "other_related_references(genome_id,species_name,radius,ANI,AF)" TEXT,
  "msa_percent" REAL,
  "translation_table" TEXT,
  "red_value" REAL,
  "warnings" TEXT,
  FOREIGN KEY (user_genome)
       REFERENCES genomes ("Accession") 
       ON UPDATE CASCADE 
       ON DELETE CASCADE 
);
CREATE INDEX "ix_gtdbtk_genome" ON "gtdbtk" ("user_genome");


CREATE TABLE "checkm" (
"Bin Id" TEXT,
  "Marker lineage" TEXT,
  "# genomes" INTEGER,
  "# markers" INTEGER,
  "# marker sets" INTEGER,
  "Completeness" REAL,
  "Contamination" REAL,
  "Strain heterogeneity" REAL,
  "Genome size (bp)" INTEGER,
  "# ambiguous bases" INTEGER,
  "# scaffolds" INTEGER,
  "# contigs" INTEGER,
  "N50 (scaffolds)" INTEGER,
  "N50 (contigs)" INTEGER,
  "Mean scaffold length (bp)" REAL,
  "Mean contig length (bp)" REAL,
  "Longest scaffold (bp)" INTEGER,
  "Longest contig (bp)" INTEGER,
  "GC" REAL,
  "GC std (scaffolds > 1kbp)" REAL,
  "Coding density" REAL,
  "Translation table" REAL,
  "# predicted genes" REAL,
  "0" TEXT,
  "1" TEXT,
  "2" TEXT,
  "3" TEXT,
  "4" TEXT,
  "5+" TEXT,
  FOREIGN KEY ("Bin Id")
       REFERENCES genomes ("Accession") 
       ON UPDATE CASCADE 
       ON DELETE CASCADE 
);
CREATE INDEX "ix_checkm_genome" ON "checkm" ("Bin Id");



