# COMEBin
GitHub repository for the manuscript "Effective binning of metagenomic contigs using COntrastive Multi-viEw representation learning".
- [Overview](#overview)
- [System Requirements](#requirements)
- [Install COMEBin via bioconda](#install)
- [Install COMEBin via source code](#started)
- [A test dataset to demo COMEBin](#demo)
- [Preprocessing](#preprocessing)
- [How to run COMEBin](#runcomebin)
- [References](#References)
- [Contacts and bug reports](#contact)
  
## <a name="overview"></a>Overview
The framework of COMEBin is shown in the following figure, which is mainly divided into the following steps: 1) Data augmentation: construct five sets of augmented data by randomly extracting subsequences from the original contigs, resulting in six views for each original contig; 2) Construct feature vector: construct nucleotide frequency feature and coverage feature for each contig (including the original sequences and augmented sequences); 3) Contrastive learning: obtain low-dimensional embedding representations suitable for binning with heterogeneous information based on multi-view contrastive learning, and 4) Clustering: generate binning results based on community division algorithm Leiden.

Among them, the network structure used in contrastive learning includes two parts: 1) "Coverage network": process coverage features and 2) "Combine network": integrate the k-mer features and the "Coverage network" to obtain representations containing heterogeneous information.

<p align="center">
<img src="https://github.com/ziyewang/COMEBin/blob/master/overview.png" width="550"/>
</p>

## <a name="requirements"></a>System Requirements
### Hardware requirements
COMEBin requires only a standard computer with enough RAM to support the in-memory operations.

### OS Requirements
COMEBin v1.1.0 is supported and tested in Linux systems.

### Dependencies
Besides the Python packages installed by `comebin_env.yaml` (which intentionally does not install PyTorch, see the install sections above), the full pipeline relies on these command-line tools:

| Tool | Where it is used |
|------|------------------|
| `bedtools` | computing per-contig coverage from each bam file (`bedtools genomecov -bga -ibam`) |
| `bowtie2` (or `bwa`) and `samtools` | preprocessing: aligning reads to contigs and producing coordinate-sorted bam files (see the preprocessing section) |
| `FragGeneScan` and `run_FragGeneScan.pl` | predicting seed genes from the contigs for marker-based clustering |
| `hmmsearch` (HMMER) | searching seed proteins against the bundled marker HMMs (`auxiliary/bacar_marker.hmm`) |
| `perl` | running the bundled marker wrapper script (`auxiliary/test_getmarker_2quarter.pl`) |
| `checkm` (CheckM v1, incl. its data bundle) | profiling candidate bins during the final selection (the `scripts/checkm_ms/*.ms` marker sets are bundled) |

All of them are available from conda-forge/bioconda, e.g.:

```sh
conda install -c conda-forge -c bioconda bedtools samtools bowtie2 bwa hmmer fraggenescan checkm-genome
```

`checkm` needs its reference data directory; on first use point it to the data, e.g. `checkm data setRoot /path/to/checkm_data`.

## <a name="install"></a>Install COMEBin via bioconda
COMEBin can be installed from Bioconda with the GPU-enabled PyTorch build from conda-forge:
```sh
conda create -n comebin_env comebin pytorch-gpu cuda-version=13.0 \
  --channel conda-forge --channel bioconda --strict-channel-priority
conda activate comebin_env
```
For a CPU-only environment, use:

```sh
conda create -n comebin_env comebin pytorch-cpu \
  --channel conda-forge --channel bioconda --strict-channel-priority
conda activate comebin_env
```

The Bioconda package declares PyTorch as a runtime dependency. The `pytorch-gpu` and `pytorch-cpu` metapackages select the corresponding conda-forge build. The Bioconda installation uses the conda-forge PyTorch 2.13 build with CUDA 13.0 by default. To use CUDA 12.6 or another CUDA version not provided by conda-forge, please install COMEBin from source and then install the corresponding CUDA-enabled PyTorch build with pip.

## <a name="started"></a>Install COMEBin via source code
You can also install COMEBin from the source code. 
After installing Anaconda (or miniconda), first, obtain COMEBin:

```sh
git clone https://github.com/ziyewang/COMEBin.git
```
Then, create an environment to run COMEBin.

```sh
cd path_to_COMEBin
conda env create -f comebin_env.yaml
conda activate COMEBin
python -m pip install torch==2.13.0 --index-url https://download.pytorch.org/whl/cu126
```

`comebin_env.yaml` installs COMEBin's other dependencies but intentionally does not install PyTorch. The command above installs PyTorch 2.13.0 with CUDA 12.6 by default. For the CUDA 13.0 build or a CPU-only build, use the [official PyTorch installation selector](https://pytorch.org/get-started/locally/).

## <a name="demo"></a>A test dataset to demo COMEBin
We provide a small dataset to demo and test the software. Test data is available at https://drive.google.com/file/d/1xWpN2z8JTaAzWW4TcOl0Lr4Y_x--Fs5s/view?usp=sharing.
The inputs for COMEBin include contigs and BAM files (reads mapping to the contigs).
```sh
Contig file: comebin_test_data/BATS_SAMN07137077_METAG.scaffolds.min500.fasta.f1k.fasta
BAM files: comebin_test_data/bamfiles/SRR5720343.bam
```
Run COMEBin on the test dataset:
```sh
CUDA_VISIBLE_DEVICES=0 bash path_to_COMEBin/COMEBin/run_comebin.sh -a path_to_comebin_test_data/BATS_SAMN07137077_METAG.scaffolds.min500.fasta.f1k.fasta \
-p path_to_comebin_test_data/bamfiles \
-o path_to_comebin_test_data/run_comebin_test \
-n 6 \
-d cuda \
-t 40
```
Excepted output is given in path_to_comebin_test_data/run_comebin_test

```sh
Final result (bins): path_to_comebin_test_data/run_comebin_test/comebin_res/comebin_res_bins
Final result in tsv format: path_to_comebin_test_data/run_comebin_test/comebin_res/comebin_res.tsv
```


## <a name="preprocessing"></a>Preprocessing

The preprocessing steps aim to generate bam files as input to our program.

Several binning methods can generate bam files by aligning reads to contigs (such as MetaWRAP), and we provide one way to generate the input files as follows.
### Generate bam files
To generate bam files from sequencing reads directly, run the script slightly modified from the "binning.sh" of MetaWRAP. The script supports different types of sequencing reads, and the default type is "paired" ([readsX_1.fastq readsX_2.fastq ...]).

```sh
cd path_to_COMEBin
cd COMEBin/scripts

bash gen_cov_file.sh -a contig_file \
-o output_dir_of_bamfiles \
path_to_sequencing_reads/*fastq

Options:

        -a STR          metagenomic assembly file
        -o STR          output directory (to save the coverage files)
        -b STR          directory for the bam files (optional)
        -t INT          number of threads (default=1)
        -m INT          amount of RAM available (default=4)
        -l INT          minimum contig length to the bin (default=1000bp).
        --single-end    non-paired reads mode (provide *.fastq files)
        --interleaved   input read files contain interleaved paired-end reads
        -f              Forward read suffix for paired reads (default="_1.fastq")
        -r              Reverse read suffix for paired reads (default="_2.fastq")

```

And the users can run the following command to keep the contigs longer than 1000bp for binning.

```sh
cd path_to_COMEBin
cd COMEBin/scripts

python Filter_tooshort.py final.contigs.fa 1000
```


## <a name="runcomebin"></a>How to run COMEBin
### Run COMEBin via bioconda
```sh
conda activate comebin_env

run_comebin.sh -a ${contig_file} \
-o ${output_path} \
-p ${path_to_bamfiles} \
-d cuda \
-t 40
```

### Run COMEBin via source code
```sh
conda activate COMEBin

bash path_to_COMEBin/COMEBin/run_comebin.sh -a ${contig_file} \
-o ${output_path} \
-p ${path_to_bamfiles} \
-d cuda \
-t 40
```

### Limiting memory and CPU usage
`run_comebin.sh` exposes several options (see `bash run_comebin.sh -h`) to fit COMEBin on machines with limited RAM or CPUs:

| Option | What it controls | Effect on resource usage |
|--------|------------------|--------------------------|
| `-d cpu` (or `cuda:N`) | training device | CPU avoids GPU memory, but the 200-epoch training is slowest on CPU |
| `-b INT` | training batch size (default 1024) | the main lever on training memory; try 384-512 on low-RAM machines |
| `-t INT` | PyTorch intra-op threads (default 5) | lower reduces the per-thread kernel workspace |
| `-n INT` | number of views (default 6) | lower reduces feature/DataLoader memory |
| `-e INT` / `-c INT` | combine/coverage network embedding sizes (default 2048) | lower reduces model parameter memory |
| `-w INT` | concurrent Leiden workers (default = `-t`) | lower reduces the peak memory of the clustering step |
| `-m INT` | HNSW neighbors kept per contig (default 100) | lower reduces index memory during clustering |
| `-s INT` | random seed for reproducible runs | -- |

Practical notes for low-memory runs:

- The training DataLoader uses a fixed worker count (currently 4) regardless of `-t`, and every worker holds a copy of the feature tensors. Prefer reducing `-b`/`-e`/`-n` over relying on thread settings alone.
- Coverage calculation (augmentation step) runs one `bedtools genomecov` per bam; it is single-threaded per bam but bam files from different samples are processed in parallel with up to `--num_threads` workers.
- Do not run several COMEBin training processes at the same time: peak memory per process does not shrink, so combined usage can easily exceed a few GiB.
- If the process is killed with status 137 (SIGKILL), it hit a memory limit (often a container/cgroup limit). Check available memory with `free -h` (and `cat /sys/fs/cgroup/memory.max` in containers) before launching.
- To additionally cap the native/OpenMP/BLAS thread pools:

```sh
export OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2
```

#### Worked example: CPU-only run inside a ~4 GiB memory cap
These numbers were measured while running the bundled test dataset (29,434 contigs, one ~5 GB bam) on a CPU-only container limited to 3.8 GiB (cgroup `memory.max`):

- **Data augmentation** finished in ~11 minutes. Coverage comes from one `bedtools genomecov` stream per bam; nothing to tune for a single-sample run.
- **Training is the bottleneck.** With `--batch_size 384 --num_threads 2` peak RSS stays under the cap:
  - main process ≈ 1.5-1.6 GB,
  - four fixed DataLoader workers ≈ 0.53 GB each,
  - combined ≈ 3.5 GB of the 3.8 GiB cap.
- **Timing** at 2 threads ≈ 8-9 s per batch, ≈ 10 min per epoch, up to ~33 h for 200 epochs. With `--earlystop` the loop breaks after top1 accuracy stays > 99% for 3 consecutive epochs (checked from epoch 10 on), which usually ends the run much earlier.

Because `main.py train` is long-running on CPU, launch it detached so it survives the launching shell being closed:

```sh
conda activate COMEBin
cd path_to_COMEBin/COMEBin
OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 PYTHONUNBUFFERED=1 \
  setsid nohup python main.py train \
    --data output_dir/data_augmentation \
    --emb_szs_forcov 2048 --emb_szs 2048 --n_views 6 \
    --batch_size 384 --num_threads 2 --earlystop --addvars --vars_sqrt \
    --device cpu --output_path output_dir/comebin_res \
    > /tmp/comebin_train.log 2>&1 < /dev/null &
```

Monitor progress with `tail -f /tmp/comebin_train.log`; the tqdm line shows the current batch, and `loss`/`top1` lines are printed every 20 steps. When it finishes, `comebin_res/` contains `embeddings.npy`, `covembeddings.npy`, `embedding_manifest.json`, etc., and the remaining stages can be run exactly as `run_comebin.sh` does:

```sh
# marker-gene seed file (auto-generated with FragGeneScan + hmmsearch if missing)
seed_file="${contig_file}.bacar_marker.2quarter_lencutoff_1001.seed"

python main.py bin --contig_file "${contig_file}" \
  --emb_file output_dir/comebin_res/embeddings.npy \
  --output_path output_dir/comebin_res \
  --seed_file "${seed_file}" --num_threads 2 --hmm_evalue 1e-5 --max_edges 100

python main.py get_result --contig_file "${contig_file}" \
  --output_path output_dir/comebin_res \
  --emb_file output_dir/comebin_res/embeddings.npy \
  --seed_file "${seed_file}" --num_threads 2 --max_edges 100 --hmm_evalue 1e-5
```

The binning result is written to `output_dir/comebin_res/comebin_res_bins` and `comebin_res.tsv` under `output_dir/comebin_res` (see the How to run section).

## <a name="References"></a>References
[1] Meyer F, Fritz A, Deng Z L, et al. Critical assessment of metagenome interpretation: the second round of challenges[J]. Nature methods, 2022, 19(4): 429-440.

[2] Parks D H, Imelfort M, Skennerton C T, et al. CheckM: assessing the quality of microbial genomes recovered from isolates, single cells, and metagenomes[J]. Genome research, 2015, 25(7): 1043-1055.

[3] https://github.com/dparks1134/UniteM.

[4] Pan S, Zhu C, Zhao X M, et al. A deep siamese neural network improves metagenome-assembled genomes in microbiome datasets across different environments[J]. Nature communications, 2022, 13(1): 2326.

## <a name="contact"></a>Contacts and bug reports
Please feel free to send bug reports or questions to
Ziye Wang: zwang17@fudan.edu.cn and Prof. Shanfeng Zhu: zhusf@fudan.edu.cn

## <a name="Citation"></a>Citation

Wang, Z., You, R., Han, H. et al. Effective binning of metagenomic contigs using contrastive multi-view representation learning. Nat Commun 15, 585 (2024). https://doi.org/10.1038/s41467-023-44290-z
