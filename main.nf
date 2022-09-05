Strategy:
1. Classify reads with Kraken2 against a comprehensive db.
  - Produce report and summary statistics.
  - Extract relevant reads - e.g. Virus
2. de novo assembly
  - Produce report.
  - Filter out "bad" contigs?
3. Taxonomic identification of contigs
  - Blast against a pre-defined database. Can also be NCBI?
  - Extract top hit from each contig.
  - Identify whether all hits are the same or not. 
4. Contiguate assembly
  - Do separately for each of the different top blast hits
