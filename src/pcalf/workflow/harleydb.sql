CREATE TABLE "gly1" ("sequence_id" TEXT,
  "sequence" TEXT
);
CREATE INDEX "ix_gly1_sequence_id" ON "gly1" ("sequence_id");

CREATE TABLE "gly2" ("sequence_id" TEXT,
  "sequence" TEXT
);
CREATE INDEX "ix_gly2_sequence_id" ON "gly2" ("sequence_id");

CREATE TABLE "gly3" ("sequence_id" TEXT,
  "sequence" TEXT
);
CREATE INDEX "ix_gly3_sequence_id" ON "gly3" ("sequence_id");

CREATE TABLE "glyx3" ("sequence_id" TEXT,
  "sequence" TEXT
);
CREATE INDEX "ix_glyx3_sequence_id" ON "glyx3" ("sequence_id");

CREATE TABLE "nterdb" ("nter" TEXT,
  "sequence_id" TEXT,
  "sequence" TEXT
);
CREATE INDEX "ix_nterdb_sequence_id" ON "nterdb" ("sequence_id");

CREATE TABLE "summary" ("sequence_id" TEXT,
  "sequence_src" TEXT,
  "flag" TEXT,
  "nter" TEXT,
  "nter_neighbor" TEXT,
  "cter" TEXT,
  "sequence" TEXT,
  "iteration" INTEGER
);
CREATE INDEX "ix_summary_sequence_id" ON "summary" ("sequence_id");

CREATE TABLE "features" ("sequence_id" TEXT,
  "sequence_src" TEXT,
  "feature_type" TEXT,
  "feature_start" INTEGER,
  "feature_end" INTEGER,
  "feature_id" TEXT,
  "pident" REAL,
  "e-value" REAL,
  "feature_src" TEXT,
  "feature_target_len" INTEGER,
  "feature_seq" TEXT
);
CREATE INDEX "ix_features_sequence_id"ON "features" ("sequence_id");

CREATE TABLE "ccya" ("sequence_id" TEXT,
  "ccyA_genomic_region" TEXT,
  "ccyA_start" INTEGER,
  "ccyA_stop" INTEGER,
  "ccyA_frame" TEXT,
  "ccyA_partial" TEXT,
  "ccyA_pseudo" TEXT,
  "ccyA_seq" TEXT
);
CREATE INDEX "ix_ccya_sequence_id"ON "ccya" ("sequence_id");

CREATE TABLE "genomes" ("Assembly Accession" TEXT,
  "latest version" TEXT,
  "Organism Name" TEXT,
  "Assembly Level" TEXT,
  "Assembly Name" TEXT,
  "is_mag" INTEGER,
  "ANI Best ANI match ANI" REAL,
  "ANI Best ANI match Assembly" TEXT,
  "ANI Best ANI match Assembly Coverage" REAL,
  "ANI Best ANI match Organism" TEXT,
  "Annotation Release Date" TEXT,
  "Assembly BioProject Lineage Accession" TEXT,
  "Assembly BioProject Lineage Title" TEXT,
  "Assembly BioSample Accession" TEXT,
  "Assembly BioSample Description Organism Taxonomic ID" INTEGER,
  "Assembly BioSample Owner Name" TEXT,
  "isolation_source" TEXT,
  "geo_loc_name" TEXT,
  "env_broad_scale" TEXT,
  "Download date" TEXT,
  "sequence_id" TEXT,
  "nter" TEXT,
  "cter" TEXT,
  "flag" TEXT,
  "is_multiple_calc" INTEGER,
  "Genome size (Mp);" TEXT,
  "# contigs" TEXT,
  "N50 (contigs);" TEXT,
  "Completeness" TEXT,
  "Contamination" TEXT,
  "Strain heterogeneity" TEXT,
  "classification" TEXT,
  "fastani_reference" TEXT
);
CREATE INDEX "ix_genomes_Assembly Accession"ON "genomes" ("Assembly Accession");