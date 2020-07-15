import os
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


def search_taxonID(table):
    return [element for element in table if
            element["GBQualifier_name"] == "db_xref" and element["GBQualifier_value"].split(":")[0] == 'taxon']


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
    # take a peek of what's in the result (ie. WebEnv, Count, etc.)
    for k in res:
        print(k, "=", res[k])
    paramEutils['WebEnv'] = res['WebEnv']  # add WebEnv and query_key to eUtils parameters to request esummary using
    paramEutils['query_key'] = res['QueryKey']  # search history (cache results) instead of using IdList
    paramEutils['rettype'] = 'xml'  # get report as xml
    paramEutils['retstart'] = 0  # get result starting at 0, top of IdList
    paramEutils['retmax'] = 5  # get next five results
    return res['IdList']

# Genbank con datos de countries sino no anda
def run(fileName, alphabet, bootstrap, aligned):
    db = 'protein' # set search to dbVar database
    if alphabet == generic_nucleotide:
        db = 'nucleotide'
    tree_phy = evolutionaryInference.fasta_to_tree(fileName, aligned, bootstrap)

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
            # output.write(record[i]['GBSeq_primary-accession'] + ": " + str(len(sourcetab))+ "\n") # this was just a check!

            # create temporary objects to collect values of the fields of interest.. note that we run a loop through 'sourcetab'.. It is probably not necessary, since there should be only one 'source' section.. but we do it this way because the 'gene' sections will be often multiple..
            temp_countries = []
            temp_lat_lons = []
            if len(sourcetab) > 0:
                for j in range(len(sourcetab)):
                    # store country information
                    if len(search_features('country', 'GBQualifier_name', sourcetab[j]['GBFeature_quals'])) > 0:
                        temp_countries.append(
                            search_features('country', 'GBQualifier_name', sourcetab[j]['GBFeature_quals'])[0][
                                'GBQualifier_value'])  # we assume that there is only one field 'organism' per each source field..
                    else:
                        temp_countries.append('NA')
                    # store latitude and longitude
                    if len(search_features('lat_lon', 'GBQualifier_name', sourcetab[j]['GBFeature_quals'])) > 0:
                        temp_lat_lons.append(
                            search_features('lat_lon', 'GBQualifier_name', sourcetab[j]['GBFeature_quals'])[0][
                                'GBQualifier_value'])  # we assume that there is only one field 'lat_lon' per each source field..
                    else:
                        temp_lat_lons.append('NA')
                countries.append(','.join(temp_countries))
                lat_lons.append(','.join(temp_lat_lons))
            else:
                countries.append('NA')
                lat_lons.append('NA')
            print('country')
            print(countries)
        countriesByGenbank[seq_record.id] = countries[0]
    draw(countriesByGenbank, tree_phy)


def draw(countriesByGenbank, tree_phy):
    geo = Nominatim(user_agent='BioLocation', timeout=2)
    plt.figure(figsize=(16, 12))
    myMap = Basemap(projection='robin', lon_0=0, lat_0=0)
    labels = []
    location = []
    for (key, value) in countriesByGenbank.items():
        temp_location = ('NA', 'NA')
        if value != 'NA':
            place = geo.geocode(value)
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
                myMap.plot(xpt, ypt, marker='o', markerfacecolor=color, markersize=str(markersize-len(colorsInMap.get((xpt, ypt)))-5))
        else:
            if xpt != 'NA':
                colorsInMap[(xpt, ypt)] = [color]
                myMap.plot(xpt, ypt, marker='o', markerfacecolor=color, markersize=str(markersize))
        nstyle = NodeStyle()
        nstyle["fgcolor"] = color.split(':')[1]
        tree_phy.get_leaves_by_name(label)[0].set_style(nstyle)
    ts = TreeStyle()
    ts.show_branch_support = True
    tree_phy.render("mytree.png", tree_style=ts)
    fig = plt.figure()
    img = mpimg.imread('mytree.png')
    imgplot = plt.imshow(img)
    os.remove('mytree.png')
    os.remove('egfr-family.phy.iqtree')
    os.remove('egfr-family.phy.contree')
    os.remove('egfr-family.phy.model.gz')
    os.remove('egfr-family.phy.splits.nex')
    os.remove('egfr-family.phy.treefile')
