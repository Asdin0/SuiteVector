# Gunicorn settings for the Railway image (see Dockerfile). Locally, from the project root:
#
#   gunicorn --config gunicorn.conf.py
#
# Railway injects PORT. WEB_CONCURRENCY can be set as a service variable to change the
# number of worker processes without a rebuild.

import os

wsgi_app = "app:server"
# Same import layout as `cd web_app && python app.py`: chdir also puts web_app on sys.path.
chdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "web_app")

bind = f"0.0.0.0:{os.environ.get('PORT', '8080')}"
workers = int(os.environ.get("WEB_CONCURRENCY", 2))
worker_class = "gthread"
threads = 4

# app.py warms the equations of state on import (symbolic differentiation with sympy).
# Loading it once in the master lets every worker fork warm and share that memory, and
# keeps the slow import out of the worker boot timeout.
preload_app = True
timeout = 120

accesslog = "-"

# Gunicorn 25.1+ opens a control socket under $HOME by default; nothing here uses it.
control_socket_disable = True
