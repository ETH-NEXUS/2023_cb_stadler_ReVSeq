#!/usr/bin/env bash

PORT=5000
echo 'ENTRYPOINT: Starting entrypoint.sh...'

python manage.py collectstatic --noinput

python manage.py wait_for_db
python manage.py makemigrations
python manage.py migrate
python manage.py initadmin
python manage.py db init





if [ "$DJANGO_DEBUG" == "True" ]; then
  python manage.py runserver 0.0.0.0:${PORT}
  # python -m uvicorn ldm.asgi:application \
  #   --log-level debug \
  #   --access-log \
  #   --workers 1 \
  #   --host 0.0.0.0 \
  #   --port ${PORT} \
  #   --reload
else
  gunicorn revseq.asgi:application \
    --log-file - \
    --workers 4 \
    --worker-class uvicorn.workers.UvicornWorker \
    --timeout 300 \
    --bind 0.0.0.0:${PORT}
fi
