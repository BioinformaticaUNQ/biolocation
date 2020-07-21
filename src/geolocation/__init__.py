import os
import pickle
from pathlib import Path
from Bio import Entrez
from Bio import SeqIO
from Bio.Alphabet import generic_nucleotide
from geopy.geocoders import Nominatim
from ete3 import Tree, NodeStyle, TreeStyle

from src import evolutionaryInference
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

# home_dir = os.path.expanduser("~")
# my_module = os.path.join(home_dir, "Anaconda3\Library\share")
project_dir = Path(__file__).parent.parent.parent
my_module = os.path.join(project_dir, "resources")
os.environ['PROJ_LIB'] = my_module

from mpl_toolkits.basemap.test import Basemap


def search_features(name, key, table):
    return [element for element in table if element[key] == name]


def listIds(seq_record, db):
    variant = seq_record.id
    Entrez.email = 'grupo6@bioinformatica.com.ar'  # provide your email address
    paramEutils = {'usehistory': 'Y', 'RetMax': '1'}  # Use Entrez search history to cache results
    # download list of GIs for the species
    print('downloading sequences for species organism:' + variant)
    # generate query to Entrez eSearch
    eSearch = Entrez.esearch(db=db, term='("' + variant + '"[Accession])', **paramEutils)
    # get eSearch result as dict object
    res = Entrez.read(eSearch)
    eSearch.close()
    # take a peek of what's in the result (ie. WebEnv, Count, etc.)
    for k in res:
        print(k, "=", res[k])
    paramEutils['WebEnv'] = res['WebEnv']  # add WebEnv and query_key to eUtils parameters to request esummary using
    paramEutils['query_key'] = res['QueryKey']  # search history (cache results) instead of using IdList
    paramEutils['rettype'] = 'xml'  # get report as xml
    paramEutils['retstart'] = 0  # get result starting at 0, top of IdList
    paramEutils['retmax'] = 5  # get next five results
    return res['IdList']


def run(fileName, alphabet, bootstrap, aligned, quantitySequences):
    '''
    # If exists countriesByGenbank.pkl
    countriesByGenbank_file = open("countriesByGenbank.pkl", "rb")
    countriesByGenbank = pickle.load(countriesByGenbank_file)
    print(countriesByGenbank)
    draw(countriesByGenbank)
    countriesByGenbank_file.close()
    '''

    db = 'protein'  # set search to dbVar database
    if alphabet == generic_nucleotide:
        db = 'nucleotide'
    tree_phy = evolutionaryInference.fasta_to_tree(fileName, aligned, bootstrap, quantitySequences)

    # location
    countriesByGenbank = {}
    for seq_record in SeqIO.parse(fileName, "fasta"):
        list_of_ids = listIds(seq_record, db)
        if len(list_of_ids) > 0:
            Entrez.email = "grupo6@bioinformatica.com.ar"
            ### read entry in xml and extract features..
            # use Entrez.efetch to read data in .xml format
            handle = Entrez.efetch(db=db, id=list_of_ids, retmode="xml")
            record = Entrez.read(handle, validate=False)
            # reset lists for values to be stored
            countries = []
            lat_lons = []
            feattab = record[0]['GBSeq_feature-table']
            refs = record[0]['GBSeq_references']
            # get and name the 'source' section (s)
            sourcetab = search_features('source', 'GBFeature_key', feattab)
            temp_countries = []
            temp_lat_lons = []
            print("sourcetab")
            print(sourcetab)
            if len(sourcetab) > 0:
                for j in range(len(sourcetab)):
                    # store country information
                    if len(search_features('country', 'GBQualifier_name', sourcetab[j]['GBFeature_quals'])) > 0:
                        temp_countries.append(
                            search_features('country', 'GBQualifier_name', sourcetab[j]['GBFeature_quals'])[0][
                                'GBQualifier_value'])
                    else:
                        temp_countries.append('NA')
                    # store latitude and longitude
                    if len(search_features('lat_lon', 'GBQualifier_name', sourcetab[j]['GBFeature_quals'])) > 0:
                        temp_lat_lons.append(
                            search_features('lat_lon', 'GBQualifier_name', sourcetab[j]['GBFeature_quals'])[0][
                                'GBQualifier_value'])
                    else:
                        temp_lat_lons.append('NA')
                countries.append(','.join(temp_countries))
                lat_lons.append(','.join(temp_lat_lons))
            else:
                countries.append('NA')
                lat_lons.append('NA')
            countriesByGenbank[seq_record.id] = countries[0]
    draw(countriesByGenbank, tree_phy)
    # draw(countriesByGenbank)

def draw(countriesByGenbank, tree_phy):
    geo = Nominatim(user_agent='BioLocation', timeout=2)
    plt.figure(figsize=(13, 12))
    myMap = Basemap(projection='robin', lon_0=0, lat_0=0)
    labels = []
    location = []
    for (key, value) in countriesByGenbank.items():
        temp_location = ('NA', 'NA')
        if value != 'NA':
            place = geo.geocode(value.split(':')[0])
            x, y = myMap(place.longitude, place.latitude)
            temp_location = (x, y)
        labels.append(key)
        location.append(temp_location)
    myMap.drawcoastlines()
    myMap.drawcountries()
    myMap.fillcontinents(color='white')
    myMap.drawmapboundary(fill_color='aqua')

    longitudes = []
    latitudes = []
    for longitude, latitude in location:
        longitudes.append(longitude)
        latitudes.append(latitude)
    colors = list(mcolors.TABLEAU_COLORS)
    colorsInMap = {}
    markersize = 15
    for label, xpt, ypt, color in zip(labels, iter(longitudes), iter(latitudes), colors):
        if colorsInMap.get((xpt, ypt)):
            colorsInMap.get((xpt, ypt)).append(color)
            if xpt != 'NA':
                myMap.plot(xpt, ypt, marker='o', markerfacecolor=color,
                           markersize=str(markersize - len(colorsInMap.get((xpt, ypt)))-1 ))
        else:
            if xpt != 'NA':
                print(xpt)
                colorsInMap[(xpt, ypt)] = [color]
                myMap.plot(xpt, ypt, marker='o', markerfacecolor=color, markersize=str(markersize))
        nstyle = NodeStyle()
        nstyle["fgcolor"] = color.split(':')[1]
        tree_phy.get_leaves_by_name(label)[0].set_style(nstyle)
    plt.annotate('Nodos no ancestrales', xy=(0, 0), xycoords='axes fraction')
    ts = TreeStyle()
    ts.show_branch_support = True
    tree_phy.get_tree_root().render("mytree.png", tree_style=ts)
    fig = plt.figure()
    img = mpimg.imread('mytree.png')
    imgplot = plt.imshow(img)
    os.remove('egfr-family.phy.iqtree')
    os.remove('egfr-family.phy.contree')
    os.remove('egfr-family.phy.model.gz')
    os.remove('egfr-family.phy.splits.nex')
    os.remove("egfr-family.phy.bionj")
    os.remove("egfr-family.phy.ckp.gz")
    os.remove("egfr-family.phy.mldist")

'''
def draw(countriesByGenbank): # tree_phy
    # countriesByGenbank_file = open("countriesByGenbank.pkl", "wb")
    # pickle.dump(countriesByGenbank, countriesByGenbank_file)
    # countriesByGenbank_file.close()
    geo = Nominatim(user_agent='BioLocation', timeout=2)
    plt.figure(figsize=(13, 12))
    myMap = Basemap(projection='robin', lon_0=0, lat_0=0)
    labels = []
    location = []
    for (key, value) in countriesByGenbank.items():
        temp_location = ('NA', 'NA')
        if value != 'NA':
            place = geo.geocode(value.split(':')[0])
            x, y = myMap(place.longitude, place.latitude)
            temp_location = (x, y)
        labels.append(key)
        location.append(temp_location)
    myMap.drawcoastlines()
    myMap.drawcountries()
    myMap.fillcontinents(color='white')
    myMap.drawmapboundary(fill_color='aqua')

    longitudes = []
    latitudes = []
    for longitude, latitude in location:
        longitudes.append(longitude)
        latitudes.append(latitude)
    colors = list(mcolors.TABLEAU_COLORS)
    colorsInMap = {}
    markersize = 15
    for label, xpt, ypt, color in zip(labels, iter(longitudes), iter(latitudes), colors):
        if colorsInMap.get((xpt, ypt)):
            colorsInMap.get((xpt, ypt)).append(color)
            if xpt != 'NA':
                myMap.plot(xpt, ypt, marker='o', markerfacecolor=color,
                           markersize=str(markersize )) # - len(colorsInMap.get((xpt, ypt)))-1
        else:
            if xpt != 'NA':
                colorsInMap[(xpt, ypt)] = [color]
                myMap.plot(xpt, ypt, marker='o', markerfacecolor=color, markersize=str(markersize))'''