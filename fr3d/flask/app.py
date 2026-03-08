from flask import Flask
from flask import jsonify
from flask import redirect
from flask import render_template
from flask import request
from flask import Response
from flask import send_from_directory
from flask import url_for
from flask import send_file

import gzip
import json
import logging
import os
import re
import subprocess
import sys
import threading
import time

logging.info("Starting BGSU-specific imports")

# from email_service import send_email_to_user

# logging for regular users ... seems to be overridden by wsgi.py
# logging.basicConfig(filename='/var/www/fr3d/app/flask.log', level=logging.DEBUG, format='%(asctime)s %(message)s')


logging.info("Setting up flask app")

app = Flask(__name__)

# turn on more detailed debugging
app.debug = True

# pipe any print statements to files
sys.stdout = open('/var/www/fr3d/flask_stdout.log', 'w')
sys.stderr = open('/var/www/fr3d/flask_stderr.log', 'w')

results_path = '/var/www/fr3d/app/results'


def create_new_filename(base_path,length=14,prefix=None):
    import random
    import string

    if not prefix:
        # get current year and month
        year = time.strftime("%Y")
        month = time.strftime("%m")
        prefix = "fr3d_" + year[2:] + month

    random_name = prefix + ''.join(random.choices(string.ascii_lowercase + string.digits, k=length-2))
    while os.path.exists(os.path.join(base_path, random_name+".json")):
        random_name = prefix + ''.join(random.choices(string.ascii_lowercase + string.digits, k=length-2))

    return random_name


def clean_query_parameters(qp):
    """
    Remove unwanted and potentially harmful characters from the input
    """

    # convert the result of request.args to a regular dictionary
    qp = qp.to_dict(flat=True)

    for key in list(qp.keys()):
        value = qp[key]
        if key == 'description':
            # allow \n
            value = value.replace("#","")     # remove explicitly
            value = value.replace("\\n","#")  # temporarily rename
            value = re.sub(r"[^A-Za-z0-9<>\|_'=~&#+!?\[\]\(\):.,\- ]", '', value)
            value = value.replace("#","\\n")  # restore \n now that all other \ are gone
        elif key == 'hs':
            value = re.sub(r'[^0-9.]', '', value)
        else:
            value = re.sub(r"[^A-Za-z0-9<>\|_'=~&+!?\[\]\(\):.,\- ]", '', value)
        qp[key] = value

    return qp


@app.route('/')
def home():
    logging.info("Trying to run the home route")
    return render_template("webfr3d.html", input_parameters={})


@app.route('/search', methods=['POST'])
def fr3d_search():
    # Get the JSON data from the request
    Q = request.get_json()

    new_name = create_new_filename(results_path,14)

    if Q.get("description","").startswith("Example ") and "}}" in Q["description"] and os.path.exists('/var/www/fr3d/app/results/making_examples.txt'):
        fields = Q["description"].split("}}")    # kind of like a password
        if len(fields) > 1:
            new_name = fields[0].lower().replace(" ","_")
            Q["description"] = Q["description"].replace("}}",":")

    # save the query data to a json file in the output directory
    OUTPUTPATH = "/var/www/fr3d/app/results/"
    JSONFILENAME = os.path.join(OUTPUTPATH, new_name+".json")
    with open(JSONFILENAME, 'w') as f:
        json.dump(Q, f, indent=2)

    logging.info("Running FR3D search with query: %s" % Q)

    # set parameters for running the FR3D search on the server
    Q["DATAPATHUNITS"] = "/var/www/html/units"       # location of NA_unit.pickle files
    Q["DATAPATHPAIRS"] = "/var/www/html/pairs"       # location of NA_pairs.pickle files
    Q["PDBDATAFILEPATH"] = "/var/www/html"           # location of NA_datafile.pickle

    # ln -s /usr/local/pipeline/fr3d-python/fr3d/search/template.html template.html
    Q["TEMPLATEFILENAME"]  = "/var/www/fr3d/app/templates/template.html"

    Q["OUTPUTPATH"]    = OUTPUTPATH
    # Q["HTMLPATH"]       = Q["OUTPUTPATH"]
    Q["JSONFILENAME"]   = JSONFILENAME
    Q["HTMLFILENAME"]   = os.path.join(Q["OUTPUTPATH"], new_name+".html")
    # Q["CSVFILENAME"]    = os.path.join(Q["OUTPUTPATH"], new_name+".csv")
    Q['url'] = "https://rna.bgsu.edu/fr3d/results/%s.html" % new_name
    Q["seeModifyQuery"] = '<a href="https://rna.bgsu.edu/fr3d/modify?id=%s">See and modify query</a> ' % new_name
    Q["gzip"] = True

    Q["MAXTIME"] = 1200             # seconds
    Q["MAXCANDIDATESHEATMAP"] = 300
    Q["MAXCANDIDATES"] = 1000
    Q["REFRESHTIME"] = 2

    # for now, since .gz files are not being served, just don't compress
    Q["gzip"] = False

    Q["PDBONLY"] = True
    Q['JSLOCATION'] = 'https://rna.bgsu.edu/rna3dhub/'

    Q['stdout'] = os.path.join(OUTPUTPATH, new_name+".stdout")
    Q['printCumulativeCandidates'] = True
    Q['downloadDataFiles'] = False           # on the server, if they don't exist, can't download them

    # write an initial output page, which will refresh every 2 seconds
    initial_page = """
    <html>
    <head>
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <title>WebFR3D search progress</title>
        <meta http-equiv="refresh" content="2">  <!-- Automatically refresh every 2 seconds -->
        <link rel="shortcut icon" href="https://rna.bgsu.edu/rna3dhub/icons/F_icon.png">
    </head>
    <body>
        <h2>FR3D output</h2>
        <p>%s</p>
        <p>This page will update every 2 seconds as results become available.</p>
        <p>If it seems to be taking too long, try a hard refresh of the page with Ctrl-F5 or similar key combination.</p>
    </body>
    </html>
    """

    if Q.get('description',''):
        description_text = "Query description: %s" % Q['description']
    else:
        description_text = "Query %s" % new_name

    with open(Q["HTMLFILENAME"], 'w') as f:
        f.write(initial_page % (description_text))

    # do a FR3D search in the background using the json query
    from fr3d.search.FR3D import fr3d_search_stdout
    task_thread = threading.Thread(target=fr3d_search_stdout, args=(Q,))
    task_thread.start()

    """
    # email the user if requested
    if Q.get("email",""):
        # basic cleanup of email addresses
        Q["email"] = Q["email"].split(",")
        if len(Q["email"]) > 5:
            # at most 5 email addresses to avoid becoming a spam server
            Q["email"] = Q["email"][:5]
        # send email to the user now, so they can look up results later
        Q = send_email_to_user(Q)
    """

    # meanwhile, direct the user to the initial page for now
    relative_html_path = os.path.join("results", new_name+".html")

    return jsonify({"status": "success", "redirect": relative_html_path})

@app.route('/results/<path:filename>')
def serve_public(filename):
    return send_from_directory(results_path, filename)


@app.route('/modify')
def modify_query():
    # example:
    query_parameters = clean_query_parameters(request.args)

    id = query_parameters.get('id','')
    if not id.endswith(".json"):
        id += ".json"
    json_file = os.path.join(results_path, id)

    if not os.path.exists(json_file):
        return "No query found with id %s" % id

    with open(json_file, 'r') as f:
        query = json.load(f)

    if "example_" in id:
        query["exampleDescription"] = query.get("description","")
        query["description"] = ""
        query["resultHTML"] = "results/%s.html" % id.replace(".json","")

    return render_template("webfr3d.html",input_parameters=query)


@app.route('/seemodify')
def see_modify_query():
    # Applies to WebFR3D results from before September 20, 2024
    # example:  https://rna.bgsu.edu/webfr3d/Results/65e00d05be4fb/65e00d05be4fb.html

    query_parameters = clean_query_parameters(request.args)

    id = query_parameters.get('id','')

    json_file = os.path.join("/var/www/html/webfr3d/Results/", id, "Query_" + id + ".json")

    if not os.path.exists(json_file):
        return "No query found with id %s" % id

    if os.path.getsize(json_file) == 0:
        return "Query file for id %s is empty" % id

    with open(json_file, 'r') as f:
        query = json.load(f)

    # map previous query parameters to the new ones
    Q = {}
    if "numpositions" in query:
        Q["numPositions"] = int(query["numpositions"])
    if "unitID" in query:
        Q["unitID"] = query["unitID"]
    if "discrepancy" in query:
        Q["discrepancy"] = query["discrepancy"]
    if "interactionMatrix" in query:
        for i in range(len(query["interactionMatrix"])):
            for j in range(len(query["interactionMatrix"][i])):
                key = 'constraintMatrix,%d,%d' % (i+1,j+1)       # count starting from 1
                if query["interactionMatrix"][i][j]:             # count starting from 0
                    Q[key] = query["interactionMatrix"][i][j]
    if "structuresToSearch" in query:
        Q["searchFiles"] = query["structuresToSearch"]
    if "repSetRelease" in query:
        Q["repSetRelease"] = query["repSetRelease"]
    if "repSetResolution" in query:
        Q["repSetResolution"] = query["repSetResolution"]
    if "name" in query and query["name"]:
        Q["description"] = query["name"]
    if "email" in query and query["email"]:
        Q["email"] = query["email"]

    return render_template("webfr3d.html",input_parameters=Q)


@app.route("/svg_viewer")
def svg_viewer():
    #
    return send_file("templates/svg_viewer.html")


@app.route('/r3dcid')
def r3dcid():
    """
    Examples
    https://rna.bgsu.edu/fr3d/r3dcid?chains=5J7L|1|AA
    """

    query_parameters = clean_query_parameters(request.args)

    logging.info("First r3dcid circular diagram request: %s" % str(query_parameters))

    if query_parameters.get('input_form','').lower() == 'true':
        return render_template("r3dcid.html",input_parameters=query_parameters)

    if query_parameters.get('chains','') == '':
        return render_template("r3dcid.html",input_parameters=query_parameters)

    if query_parameters.get('format','') == 'html':
        return render_template("svg_viewer.html",input_parameters=query_parameters)
        return redirect(url_for("svg_viewer", **query_parameters))

    from flask import send_file
    import os
    import r3dcid

    chains_string = query_parameters.get('chains','')

    display_format = query_parameters.get('format','pdf').lower()
    if not display_format in ['pdf','svg']:
        display_format = 'pdf'

    # set parameters if they are specified in the URL
    params = {}
    if 'coloring' in query_parameters:
        params['coloring'] = query_parameters.get('coloring','default')
    if 'show' in query_parameters:
        params['show'] = query_parameters.get('show','')
    if 'dim' in query_parameters:
        params['dim'] = query_parameters.get('dim','')
    if 'hide' in query_parameters:
        params['hide'] = query_parameters.get('hide','')
    if 'text' in query_parameters:
        params['text'] = query_parameters.get('text','basepair')
    if 'n3d' in query_parameters:
        params['n3d'] = query_parameters.get('n3d',True)
        if params['n3d'].lower() == 'false':
            params['n3d'] = False
        elif params['n3d'].lower() == 'true':
            params['n3d'] = True
    if 'header' in query_parameters:
        params['header'] = query_parameters.get('header','')
    if 'description' in query_parameters:
        d = query_parameters.get('description','')
        if len(d) > 300:
            d = d[0:299] + "+"
        if len(d) > 0:
            params['description'] = d

    if 'assemblies' in query_parameters:
        params['assembly'] = query_parameters.get('assemblies','')

    if 'symmetries' in query_parameters:
        params['symmetry'] = query_parameters.get('symmetries','')

    hs = query_parameters.get('hs','-1')
    try:
        hs = int(hs)
    except:
        hs = -1
    if hs >= 0:
        params['helix_size'] = hs

    output_path = '/var/www/fr3d/app/r3dcid_output'
    params['data_directory'] = '/var/www/html/pairs/'
    params['output_directory'] = output_path

    logging.info("params %s" % params)

    try:
        if 'description' in params:
            # generate PDF or SVG output, not both, with this description
            params['format'] = display_format
            # get the filename without extension
            filename, message = r3dcid.main(chains_string, params)
        else:
            # create both pdf and svg formatted output with no description
            params['format'] = 'pdf,svg'

            # check to see if the file already exists, to avoid generating it again
            filename = r3dcid.get_filename(chains_string, params)

            pdf_file = os.path.join(output_path,filename+".pdf")
            if not os.path.exists(pdf_file):
                filename, message = r3dcid.main(chains_string, params)
            else:
                svg_gz_file = os.path.join(output_path,filename+".svg.gz")
                if not os.path.exists(svg_gz_file):
                    filename, message = r3dcid.main(chains_string, params)

    except Exception as e:
        exception_type, exception_object, exception_traceback = sys.exc_info()
        # line_number = exception_traceback.tb_lineno
        output = ""
        output += "\nSomething went wrong with this request: %s\n" % (exception_type)
        output += "%s\n" % type(e)
        output += "%s\n" % exception_traceback
        output += "%s\n" % e
        return output

    if len(filename) == 0:
        return "No chains found for the requested structure " + message

    # if a .ps file was produced, make sure the .pdf file is made
    pdf_file = os.path.join(output_path,filename+".pdf")
    ps_file  = os.path.join(output_path,filename+".ps")
    if os.path.exists(ps_file) and not os.path.exists(pdf_file):
        try:
            subprocess.run(["ps2pdf", "-dFIXEDMEDIA", output_path, ps_file],check=True)
            os.remove(ps_file)
        except:
            raise RuntimeError("ps2pdf failed for %s" % ps_file) from e

    # gzip every .svg file
    svg_file    = os.path.join(output_path,filename+".svg")
    svg_gz_file = os.path.join(output_path,filename+".svg.gz")
    if os.path.exists(svg_file):
        try:
            subprocess.run(["gzip", "-f", svg_file],check=True)
        except subprocess.CalledProcessError as e:
            raise RuntimeError("gzip failed for %s" % svg_file) from e

    if display_format == 'pdf':
        try:
            return send_file(pdf_file, attachment_filename=filename+".pdf", as_attachment=True)  # instant download
        except Exception as e:
            return str(e)
    else:
        try:
            with gzip.open(svg_gz_file, "rb") as f:
                svg_bytes = f.read()
            return Response(svg_bytes, mimetype="image/svg+xml", headers={"Content-Disposition": f'inline; filename="{filename}.svg"'})

            # return send_file(svg_file, mimetype="image/svg+xml", download_name=filename+".svg", as_attachment=False)
        except Exception as e:
            return str(e)


# Handle errors of different types
@app.errorhandler(400)
def not_found_error(error):
    return render_template('error.html'),400
@app.errorhandler(403)
def error_403(error):
    return render_template('error.html'),403
@app.errorhandler(404)
def error_404(error):
    return render_template('error.html'),404
@app.errorhandler(429)
def error_429(error):
    return render_template('error.html'),429
@app.errorhandler(500)
def internal_error(error):
    return render_template('error.html'),500
@app.errorhandler(502)
def error_502(error):
    return render_template('error.html'),502
@app.errorhandler(503)
def error_503(error):
    return render_template('error.html'),503
@app.errorhandler(504)
def error_504(error):
    return render_template('error.html'),504


@app.route('/circular')
def redirect_to_r3dcid():
    return redirect(url_for('r3dcid'))


if __name__ == '__main__':
    # app.debug = False
    app.run()