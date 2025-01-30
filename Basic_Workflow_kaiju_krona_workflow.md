## Basic Workflow: Running Kaiju and KronaTools

**Step 1: Run Kaiju for Taxonomic Classification**

```bash
kaiju -t nodes.dmp -f kaiju_db.fmi -i reads.fastq -o kaiju_output.txt -z 8
```
- `-t nodes.dmp`: Taxonomy nodes file
- `-f kaiju_db.fmi`: Kaiju database file
- `-i reads.fastq`: Input metagenomic reads
- `-o kaiju_output.txt`: Output file
- `-z 8`: Number of CPU threads

**Step 2: Convert Kaiju Output to Krona Format**

```bash
kaiju2krona -t nodes.dmp -n names.dmp -i kaiju_output.txt -o kaiju_krona.txt
```

**Step 3: Generate Krona HTML Visualization**

```bash
ktImportText -o kaiju_krona.html kaiju_krona.txt
```

**Step 4: Open the Interactive Krona Plot**

```bash
firefox kaiju_krona.html
```
*Alternatively, open in Google Chrome or any web browser.*
