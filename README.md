# CADDIE
This app consists of a Django REST API backend with a postgresql database, a Redis Task Task Queue server and an Angular frontend.
The app is dockerized, which means the Django backend, the postgresql database and the Redis task server run in one container 
with their own images. The frontend runs in a separate container. The docker configuration for the backend 
can be found in the Dockerfile and the .gitlab-ci.yml. 


## Local deployment

Make sure that the ports in 'docker-compose.yml' are correct and it is recommend to adjust the default passwords in 'docker-django.env'.  

`docker-compose up --build -d`



# Development

### Setup

Recommended are at least 16GB RAM, everything below could cause problems while filling up the database,
especially on calculating the minimal distances since the networks to calculate them on are quite big.
Populating the database can take ~2 days.

First, setup the python environment (any python 3.6 - 3.8 will do):

`python -m venv .env`

Enter the environment via:

`.env/Scripts/activate` (windows)

`source .env/bin/activate` (linux & mac)

Install all necessary packages

`pip3 install -r requirements.txt`

---

The database is filled based on source files and pre-computed files in the data directory 
The datbase is filled automatically on docker build 
or by triggering the file "import-data.sh" (or "import-data.bat"). 
In case changes have been made to the files in the data directory, update the pre-computed files.
Note that this is only necessary if there were changes in the database or the files
compared to how it is in the main branch. Otherwise the pre-computed
files in the repo are up to date.

Networks are pre-computed for each combination of cancer driver gene dataset, gene-gene interaction
 dataset and drug-gene interaction dataset. To update the graph-tool files run:
`python3 manage.py make_graphs`



The nodes for the minimum spanning trees are pre-computed for each cancer type based on information from the database. 
To update the minimum spanning trees:

`python3 manage.py populate_db --clear_model MinSpanningTree`

`python3 manage.py make_min_span_trees`



The shortest distances are pre-computed as well. To update the shortest distances, run 

`python3 manage.py make_shortest_distance_files`

Afterwards you will need to update the database on the shortest distances with

`python3 manage.py populate_db --clear_model ShortestDistanceGeneToCancerGene,ShortestDistanceDrugToCancerGene`

`python3 manage.py populate_db --add_data shortest_distances`

The django server will run on the localhost port 8000.


### Code Structure

The code structure follows a basic django app (https://www.djangoproject.com/). Django is used in this project as a 
REST API, providing endpoints to the Angular frontend. 

All data resources are stored in the data directory. Where possible, the original data files are stored without editing.
This allows to simply exchange the source file to update CADDIE. Of course, the database will have to be updated too. 
For some sources, the files required processing. Jupyter notebook scripts which where used for the processing can be found in the 
"caddie/scratch" repository. 

Functions to load the information into the database are stored in "caddie/management/includes". 

The DataCleaner is a class
with functions to load data for each single resource. The data is processed here as well and then returned in a format,
which allows to insert the data directly into the database.

The DatabaseController provides the functionality to insert all kinds of data into the database. It receives processed data from 
the DatabaseCleaner and fills the database. 

The DatabasePopulator manages the loading and filling of the database. It calls functions from DataCleaner and 
DatabaseController.

The endpoints are defined in "caddie/views.py". All requests are handled here, received data is processed and stored 
and, upon request, data is also loaded, processed and send to the frontend.

The internal task server runs with Redis Task Queue (https://redislabs.com/ebook/part-2-core-concepts/chapter-6-application-components-in-redis/6-4-task-queues/). 
When the user starts an analysis task, the django backend receives the information, creates a task objects and runs the task 
on the Redis server. The code for the task execution can be found in the "task" directory. When a task is finished, the results
are stored in the database. The frontend sends periodically requests to the backend to check for the task status. As soon as Django 
finds the task results in the database, it loads and processes them in view.py and sends them to the frontend.


The network calculations run mostly with graph tools. Useful information about graph tool that might help debugging can be found here:

`https://git.skewed.de/count0/graph-tool/-/wikis/installation-instructions`



### Cleanup tasks

To check the code, enter python environment and run:

`pycodestyle`

`flake8`



### Docker PROD environment (building is optional)
``docker-compose up --build``


### Docker DEV environemt (building is optional)
``docker-compose -f docker-compose.yml -f docker-compose.dev.yml up -d --build``