## Load initial data

```bash
$ python manage.py fixture strain_substrain
$ python manage.py fixture panel_strain
$ python manage.py fixture taxon_id
$ python manage.py fixture file_type

 
```

## By new adjustments to the table

Make changes to 

API: 
1. core/models.py SampleCount
2. core/serializers.py SampleCountSerializer
3.  If it is a numeric value, the aggregate methode in core/views.py SampleCountViewSet
4. api/app/config.config.yml
5. core/management/commands/import.py counts method 

UI:
1. src/models.core.js
2. src/components/helpers.ts
3. src/components/data.ts columnsSampleCount



## To filter controls
http://localhost:8010/api/samples/?control=true

```bash

## Import



```
## Transfer data from server

To transfer only the files that do not start with "sample_" from the remote directory, you can use a combination of scp and rsync commands. 
Unfortunately, scp alone does not support excluding files based on patterns, but rsync does.

```bash
rsync -avz --exclude='sample_*' aesche@revseq.nexus.ethz.ch:/data/revseq/results/gather_results/RVSeqPlate10/ /Users/aesche/git/2023_cb_stadler_ReVSeq/data/RVSeqPlate10/
    
```

docker compose  -f docker-compose.yml -f  docker-compose.dev.yml  up
 

scp aesche@revseq.nexus.ethz.ch:/data/revseq/results/gather_results/RVSeqPlate10/RVSeqPlate10_empty_samples.txt /Users/aesche/git/2023_cb_stadler_ReVSeq/data/RVSeqPlate10/