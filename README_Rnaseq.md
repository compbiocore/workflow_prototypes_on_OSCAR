### 1. Pull Down the Repo

If you haven't already pulled down the repo

```bash 
cd $HOME
git clone https://github.com/compbiocore/workflows_on_OSCAR.git
```

### 2. Build the Docker
```bash
cd metadata
docker build -f RNASeq_Dockerfile -t rnaseq_pipeline --no-cache --platform linux/amd64 .
```