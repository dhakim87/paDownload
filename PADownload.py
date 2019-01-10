__author__ = "Daniel Hakim, Ideker Lab"

#Depends on apache requests, use "pip install requests"

import requests
import xml.etree.ElementTree as ET
import os
import sys
from urlparse import urlparse
import urllib
import uuid

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

    downloadImagesStreaming(outputXML)

#Ugh, the full protein atlas xml is 8 gigs, discard each entry after handling it.
def downloadImagesStreaming(xmlFile):
    root = None;
    nsmap = {}
    tagStack = []
    activeEntry = False;
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
                handleEntry(elem);
                activeEntry = False;
            if activeEntry == False:
                #Wow this is a huge xml file.
                elem.clear();
                root.clear();
            
    print nsmap;

def handleEntry(entry):
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
                    for assayImage in data.findall("./assayImage"):
                        for imageUrl in assayImage.findall("./image/imageUrl"):
                            print("URL: " + imageUrl.text);
                            #TODO FIXME HACK:  Figure out what metadata we want to keep with each image.  Should probably make a simple sql db or something.
                            downloadImage(imageUrl.text);


#Pull down all of the human cell line assays.  Modify this function to pull a different set of images.
def downloadImages(xmlFile):
    print("Parsing xml (This may take a few minutes for large numbers of search results)...")
    tree = ET.parse(xmlFile)
    print("Retrieving Images...");
    root = tree.getroot()
    for entry in root.findall("./entry"):
        handleEntry(entry);

def downloadImage(url):
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

for arg in sys.argv[1:]:
    ss = arg.split("=")
    if len(ss) >= 2:
        if ss[0] == "ID":
            GENE_ID = ss[1]
        elif ss[0] == "SEARCH":
            SEARCH = ss[1]
        elif ss[0] == "XML":
            XML = ss[1]

if GENE_ID is not None:
    proteinAtlasGet(GENE_ID)
elif SEARCH is not None:
    proteinAtlasSearch(SEARCH);
elif XML is not None:
    downloadImages(XML);
else:
    print "Usage: "
    print "python " + sys.argv[0] + " ID=<ensemblGeneID>"
    print "-OR-"
    print "python " + sys.argv[0] + " SEARCH=<search string>"
    print "-OR-"
    print "python " + sys.argv[0] + " XML=<pathToLocalXMLFile>"
    print ""
    print "Ex: python " + sys.argv[0] + " ID=ENSG00000134057"
    print "Ex: python " + sys.argv[0] + " SEARCH=\"prognostic:Breast cancer\""
