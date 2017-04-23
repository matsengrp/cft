
CFT Web application

## Installation

Run `pip install -r requirements.txt` (add `--user` if you get permissions errors).

Then install miniconda and ete using ete's flavor of instructions [here](http://etetoolkit.org/download/).

## Usage

```
cd cft/cftweb
python -m cftweb /path/to/data/files/*/metadata.json
# Typically this will be the following:
#python -m cftweb ../output/*/metadata.json
```

You may specify multiple json files in this fashion as long as each has a unique `dataset_id` attribute.

The default port is `5000`.
If someone else is running the web server on the same machine (or something else using that port), you can set a different one using the `-P` flag.

### Running CFT web in production

The default dev server shipped with Flask is rather stupid and can't handle multiple web requests at the same time, and also gets thrown for a loop if a socket connection closes before all the data has been sent, periodically crashing the app.
If you install `gevent` via pip or conda, the `gevent.wsgi` module is used instead, which should resolve these issues.
Eventually we'll have a more clever Docker based server setup, but for now this will keep us limping along.

Also note that for production, you may want to specify either the `--email` or `--slack` flags to the invocation above so that folks can be notified of errors.
For slack notifications (recommended), you will also need to obtain an API token (see <https://api.slack.com/web#authentication>), and set the `SLACK_TOKEN` environment variable accordingly.

