__author__ = "Daniel Hakim, Ideker Lab"

#Depends on apache requests, use "pip install requests"

import requests
import xml.etree.ElementTree as ET
import os
import sys
from urlparse import urlparse
import urllib
import uuid
import sqlite3


#Example EnsemblGeneID Gene ID: ENSG00000134057
def proteinAtlasGet(ensemblGeneID):
    url = 'http://www.proteinatlas.org/' + urllib.quote_plus(ensemblGeneID) + '.xml'
    downloadXML(url, ensemblGeneID)

#To make search strings, go to proteinatlas.org and check out their search tool.
#Example searchString: "subcell_location:Nucleoli,Nucleoli fibrillar center;Enhanced"
#Example searchString: "prognostic:Breast cancer"
def proteinAtlasSearch(searchString):
    url = 'http://www.proteinatlas.org/search/' + urllib.quote_plus(searchString) + "?format=xml&compress=no"
    downloadXML(url, uuid.uuid5(uuid.NAMESPACE_URL, url).hex)

def downloadXML(url, outputFileID):
    print("Retrieving " + url);
    outputDir = os.path.join("./output/")
    if not os.path.exists(outputDir):
        os.makedirs(outputDir)
    
    outputXML = os.path.join(outputDir, outputFileID + '.xml')
    print(" -> " + outputXML);
    #Check if xml file is already downloaded
    #TODO: add a force option
    if (not os.path.exists(outputXML)):
        # if not already downloaded, download it.
        resp = requests.get(url)
        # saving the xml file
        with open((outputXML), 'w') as f:
            f.write(resp.content)
    else:
        print("Already Downloaded, skip!");

    downloadImagesStreamingXML(outputXML)

def downloadImagesDB(dbPath, proteinListPath):
    urlsToDownload = []
    if proteinListPath is None:
        #Download all urls
        conn = sqlite3.connect(dbPath)
        cur = conn.cursor()
        cur.execute("SELECT url FROM image")
        rows = cur.fetchall()
        for row in rows:
            urlsToDownload.append(row[0])
        conn.close()
    else:
        proteinList = []
        with open(proteinListPath) as proteinFile:
            for line in proteinFile:
                proteinList.append(line.strip())
    
        print "Proteins: ", len(proteinList)
    
        #Download urls associated with named proteins
        conn = sqlite3.connect(dbPath)
        for protein in proteinList:
            cur = conn.cursor()
            cur.execute("SELECT url FROM image WHERE protein=?", (protein,))
            rows = cur.fetchall()
            for row in rows:
                urlsToDownload.append(row[0])
        conn.close()

    for url in urlsToDownload:
        downloadImage(url)

#Ugh, the full protein atlas xml is 8 gigs, discard each entry after handling it.
def downloadImagesStreamingXML(xmlFile):
    root = None;
    nsmap = {}
    tagStack = []
    activeEntry = False;
    conn = sqlite3.connect('images.db')
    c = conn.cursor()

    # Create table
    c.execute('''CREATE TABLE IF NOT EXISTS image
                 (url TEXT, protein TEXT, antibody TEXT, cell_line TEXT, location TEXT)''')
    c.execute('''CREATE UNIQUE INDEX IF NOT EXISTS image_no_dups ON image(url, protein, antibody, cell_line, location)''')
    
    try:
        for event, elem in ET.iterparse(xmlFile, events=('start', 'end', 'start-ns', 'end-ns')):
            if event == 'start-ns':
                ns, url = elem
                nsmap[ns] = url;
            if event == 'start':
    #            print("START: " + elem.tag);
                if root is None:
                    root = elem;
                tagStack.append(elem.tag);
                if elem.tag == "entry":
                    activeEntry = True;
            
            if event == 'end':
    #            print("END  : " + elem.tag);
                if tagStack[-1] != elem.tag:
                    print("BAD XML:  End Tag: " + elem.tag + " But expected: " + tagStack[-1]);
                    raise error;
                tagStack.pop();

                if elem.tag == "entry":
                    handleEntry(elem, c);
                    activeEntry = False;
                if activeEntry == False:
                    #Wow this is a huge xml file.
                    elem.clear();
                    root.clear();
    except KeyboardInterrupt:
        print "Canceled By User";
        pass;

    conn.commit()
    conn.close()

    print nsmap;

def handleEntry(entry, cursor):
    print("Name: " + entry.find("./name").text)
    db = entry.find("./identifier").get("db")
    id = entry.find("./identifier").get("id")

    if db != "Ensembl":
        print("UNSUPPORTED DB ID: " + db)
        return;

    for antibody in entry.findall("./antibody"):
        print("Antibody ID: " + antibody.get("id"))
        for cellExpression in antibody.findall("./cellExpression"):
            source = cellExpression.get("source");
            technology = cellExpression.get("technology");
            print("Source: " + source + " Technology: " + technology)
            for subAssay in cellExpression.findall("./subAssay"):
                assayType = subAssay.get("type")
                if assayType != "human":
                    print("Skipping " + assayType + " assay.")
                    continue
                for data in subAssay.findall("./data"):
                    print("Cell Line: " + data.find("cellLine").text)
                    location = data.find("location")
                    if location is None:
                        location = ""
                    else:
                        location = location.text
                    print("Location: " + location)
                    for assayImage in data.findall("./assayImage"):
                        for imageUrl in assayImage.findall("./image/imageUrl"):
                            print("URL: " + imageUrl.text);
                            downloadImage(imageUrl.text);
                            tuple = (imageUrl.text, entry.find("./name").text, antibody.get("id"), data.find("cellLine").text, location)
                            #(url TEXT, protein TEXT, antibody TEXT, cell_line TEXT)''')
                            cursor.execute("INSERT OR IGNORE INTO image VALUES (?,?,?,?,?)", tuple)

def downloadImage(url):
    if not SHOULD_DOWNLOAD:
        print("Downloads Disabled, Skipping download of: " + url)
        return

    #TODO FIXME HACK:  May need a force in case of partial image downloads.
    parsedTuple = urlparse(url);
    relativePath = parsedTuple[2];
    outputFile = "./output" + relativePath;
    print("Output File: " + outputFile);
    
    # if not already downloaded, download it.
    if (not os.path.exists(outputFile)):
        outputDir = os.path.dirname(outputFile);

        if not os.path.exists(outputDir):
            os.makedirs(outputDir)

        with open(outputFile, 'w') as f:
            resp = requests.get(url)
            f.write(resp.content)
    else:
        print("Already downloaded, skip");

GENE_ID = None
SEARCH = None
XML = None
DB = None
PROTEIN_LIST = None

SHOULD_DOWNLOAD=True

for arg in sys.argv[1:]:
    ss = arg.split("=")
    if len(ss) >= 2:
        if ss[0] == "ID":
            GENE_ID = ss[1]
        elif ss[0] == "SEARCH":
            SEARCH = ss[1]
        elif ss[0] == "XML":
            XML = ss[1]
        elif ss[0] == "DB":
            DB = ss[1]
        elif ss[0] == "PROTEIN_LIST":
            PROTEIN_LIST = ss[1]
        else:
            raise Exception("Unknown Command: " + arg);
    elif arg == "--no-download":
        SHOULD_DOWNLOAD = False
    else:
        raise Exception("Unknown Command: " + arg);

if GENE_ID is not None:
    proteinAtlasGet(GENE_ID)
elif SEARCH is not None:
    proteinAtlasSearch(SEARCH);
elif XML is not None:
    downloadImagesStreamingXML(XML);
elif DB is not None:
    downloadImagesDB(DB, PROTEIN_LIST)
else:
    print "Usage: "
    print "python " + sys.argv[0] + " ID=<ensemblGeneID>"
    print "-OR-"
    print "python " + sys.argv[0] + " SEARCH=<search string>"
    print "-OR-"
    print "python " + sys.argv[0] + " XML=<pathToLocalXMLFile>"
    print "-OR-"
    print "python " + sys.argv[0] + " DB=<pathToLocalDBFile> PROTEIN_LIST=<pathToNewlineSeparateProteinNames>"
    print ""
    print "Ex: python " + sys.argv[0] + " ID=ENSG00000134057"
    print "Ex: python " + sys.argv[0] + " SEARCH=\"prognostic:Breast cancer\""
    print ""
    print "Options: "
    print "--no-download    //This will build the images.db database but will not download any imagery"

