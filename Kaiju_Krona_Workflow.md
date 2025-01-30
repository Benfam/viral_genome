# Comprehensive Workflow: Running Kaiju and KronaTools

## 1. Install Kaiju

**Option A: Using Bioconda**

If you have Conda installed:

```bash
conda install -c bioconda kaiju
```

**Option B: Compiling from Source**

1. **Clone the Repository:**

   ```bash
   git clone https://github.com/bioinformatics-centre/kaiju.git
   cd kaiju/src
   ```

2. **Install Dependencies:**

   Ensure the zlib development library is installed. On Ubuntu:

   ```bash
   sudo apt install libz-dev
   ```

3. **Compile Kaiju:**

   ```bash
   make
   ```

   After compilation, the executables will be in the `kaiju/bin` directory. Add this directory to your system's `$PATH` or move the executables to a directory already in your `$PATH`.

## 2. Prepare the Kaiju Database

Kaiju requires a reference database for classification.

**Option A: Download Pre-built Index**

1. **Download the Index:**

   ```bash
   wget https://kaiju.binf.ku.dk/database/kaiju_db_nr_2023-05-10.tgz
   ```

2. **Extract the Files:**

   ```bash
   tar -xzf kaiju_db_nr_2023-05-10.tgz
   ```

   This will yield:

   - `kaiju_db_nr.fmi` (the index)
   - `nodes.dmp` (taxonomy nodes)
   - `names.dmp` (taxonomy names)

**Option B: Build the Index Locally**

1. **Create a Directory:**

   ```bash
   mkdir kaijudb
   cd kaijudb
   ```

2. **Build the Database:**

   ```bash
   kaiju-makedb -s refseq
   ```

   Replace `refseq` with your desired database source. Note: Building the database requires significant storage and computational resources.

## 3. Run Kaiju for Taxonomic Classification

**For Single-End Reads:**

```bash
kaiju -t nodes.dmp -f kaiju_db_nr.fmi -i reads.fastq.gz -o kaiju_output.txt -z 8
```

**For Paired-End Reads:**

```bash
kaiju -t nodes.dmp -f kaiju_db_nr.fmi -i reads_R1.fastq.gz -j reads_R2.fastq.gz -o kaiju_output.txt -z 8
```

- `-t`: Path to `nodes.dmp`
- `-f`: Path to `kaiju_db_nr.fmi`
- `-i`: Input file for single-end reads or forward reads of paired-end
- `-j`: Input file for reverse reads of paired-end
- `-o`: Output file
- `-z`: Number of CPU threads to use

## 4. Install KronaTools

To generate interactive visualizations, install KronaTools, which includes the `ktImportText` command.

**Option A: Using Bioconda**

If you have Conda installed:

```bash
conda install -c bioconda krona
```

**Option B: Compiling from Source**

1. **Clone the Repository:**

   ```bash
   git clone https://github.com/marbl/Krona.git
   cd Krona/KronaTools
   ```

2. **Run the Installation Script:**

   ```bash
   ./install.pl
   ```

   This script will install KronaTools and its components. Ensure that the installation directory is added to your system's `$PATH` to access `ktImportText` from anywhere.

**Verify Installation:**

After installation, confirm that `ktImportText` is accessible:

```bash
ktImportText -h
```

This command should display the help information for `ktImportText`, indicating a successful installation.

## 5. Convert Kaiju Output for Krona Visualization

1. **Convert Kaiju Output to Krona Format:**

   ```bash
   kaiju2krona -t nodes.dmp -n names.dmp -i kaiju_output.txt -o kaiju_krona.txt
   ```

   - `-t`: Path to `nodes.dmp`
   - `-n`: Path to `names.dmp`
   - `-i`: Kaiju output file
   - `-o`: Output file for Krona

2. **Generate Krona HTML Visualization:**

   ```bash
   ktImportText -o kaiju_krona.html kaiju_krona.txt
   ```

3. **View the Visualization:**

   Open `kaiju_krona.html` in a web browser to interactively explore your taxonomic classification results.

**Note:** Ensure all paths to files and directories are correctly specified based on your system's configuration. For comprehensive details and troubleshooting, refer to the [Kaiju GitHub repository](https://github.com/bioinformatics-centre/kaiju) and the [KronaTools GitHub repository](https://github.com/marbl/Krona).
