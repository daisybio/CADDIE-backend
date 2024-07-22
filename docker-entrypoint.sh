#!/bin/bash


python3 manage.py migrate --run-syncdb
# python3 manage.py createfixtures

python3 manage.py cleanuptasks

file="docker-entrypoint.lock"
# exit if entrypoint.lock exists to prevent new import of data every time docker is restarted
if ! test -f "$file"; then
  # create docker-entrypoint.lock to see that data was already imported
  # used for development, otherwise db reset after each docker-start at system start
  # to reload db on docker restart, simply delete the file

  if [ -z "$DB_UPDATE_ON_START" ] || [ "$DB_UPDATE_ON_START" = "0" ]
  then
    echo "Update on startup disabled!"
  else
    sh import-data.sh

    echo "building graphs"
    python3 manage.py make_graphs

    echo "building minimum spanning trees"
    python3 manage.py make_min_span_trees
    #  python manage.py make_shortest_distance_files
    #  python manage.py populate_db --add_data shortest_distances

    # add polyphen2 data for SNP lookup
    mkdir data/poly phen2/
    wget https://cloud.uni-hamburg.de/s/4dDob7rNisAjd2Q/download/polyphen-2.2.2-whess-2011_12.sqlite.bz2 -O data/polyphen2/polyphen-2.2.2-whess-2011_12.sqlite.bz2
    # unpack the bz2 sqlite database and delete compressed file
    bzip2 -d data/polyphen2/polyphen-2.2.2-whess-2011_12.sqlite.bz2
    chmod 777 data/polyphen2/polyphen-2.2.2-whess-2011_12.sqlite
  fi

  touch $file
fi

/usr/bin/supervisord -c "/etc/supervisor/conf.d/supervisord.conf"
