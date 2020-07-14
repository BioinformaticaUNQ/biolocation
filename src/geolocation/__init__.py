import os
import subprocess
from pathlib import Path
from Bio import Entrez, AlignIO, Phylo
from Bio import SeqIO
from geopy.geocoders import Nominatim
from ete3 import Tree, NodeStyle, TreeStyle
from src.evolutionaryInference import clustalo
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


def listIds(seq_record):
    variant = seq_record.id
    Entrez.email = 'grupo6@bioinformatica.com.ar'  # provide your email address
    db = 'nucleotide'  # set search to dbVar database
    # db = 'protein'
    # db = 'pubmed'
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


def dataset(fileName, alphabet, bootstrap, aligned):
    out_file = os.path.basename(fileName) + '-aligned.fasta'
    clustalo.runClustalO("grupo6@bioinformatica.com", fileName, outfilename=out_file, fmt='clustal')
    # aligned_file = open(out_file, 'r')
    # align = AlignIO.read(aligned_file, "clustal")
    # aligned_file.close()
    AlignIO.convert(out_file, "clustal", "egfr-family.phy", "phylip-relaxed")
    os.remove(out_file)
    # calculator = DistanceCalculator('blastn')
    # dm = calculator.get_distance(align)
    # constructor = DistanceTreeConstructor()
    # tree = constructor.nj(dm)
    iqtree_exe = os.path.join(project_dir, 'resources/iqtree-1.6.12-Windows/bin/iqtree.exe')
    # subprocess.run(
    #    [iqtree_exe, "-s", "egfr-family.phy", "-m", "TIM3+R5", "-alrt", "1000", "-bo", "100", "-wbtl", "-nt", "AUTO",
    #     "-redo"])
    # subprocess.call(['iqtree','-s', file , '-bb', str(bootstrap)])
    subprocess.run([iqtree_exe, '-s', "egfr-family.phy", '-bb', str(bootstrap)])
    os.remove('egfr-family.phy')
    # seqTree = open("egfr-family.phy.bionj", "r")
    seqTree = open("egfr-family.phy.treefile", "r")
    tree_phy = Tree(str(seqTree.readlines().__getitem__(0)))
    seqTree.close()
    os.remove("egfr-family.phy.bionj")
    os.remove("egfr-family.phy.ckp.gz")
    # os.remove("egfr-family.phy.boottrees")
    os.remove("egfr-family.phy.mldist")
    # os.remove("egfr-family.phy.treefile")
    # Tree.convert_to_ultrametric(tree_phy)
    # tree.rooted = True
    # tree_phy = tree.as_phyloxml()
    # tree_phy = Phylogeny.from_tree(tree_phy)
    # tree_phy.root.color = (128, 128, 128)

    # location
    countriesByGenbank = {}
    for seq_record in SeqIO.parse(fileName, "fasta"):
        list_of_ids = listIds(seq_record)
        if len(list_of_ids) > 0:
            Entrez.email = "grupo6@bioinformatica.com.ar"
            ### read entry in xml and extract features..
            # use Entrez.efetch to read data in .xml format
            handle = Entrez.efetch(db="nucleotide", id=list_of_ids, retmode="xml")
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
    # fig, ax = plt.subplots()
    myMap = Basemap(projection='robin', lon_0=0, lat_0=0)  # ax=ax
    labels = []
    location = []
    for (key, value) in countriesByGenbank.items():
        place = geo.geocode(value)
        x, y = myMap(place.longitude, place.latitude)
        labels.append(key)
        location.append((x, y))
    myMap.drawcoastlines()
    myMap.drawcountries()
    myMap.fillcontinents(color='white')
    myMap.drawmapboundary(fill_color='aqua')
    # myMap.bluemarble()
    # ax.text(x, y, 'label', ha='center', size=20)
    # myMap.plot(x, y, color='g', marker='o', markersize='15')
    # ax.scatter(x, y)

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
            myMap.plot(xpt, ypt, marker='o', markerfacecolor=color, markersize=str(markersize-len(colorsInMap.get((xpt, ypt)))-5))
        else:
            colorsInMap[(xpt, ypt)] = [color]
            myMap.plot(xpt, ypt, marker='o', markerfacecolor=color, markersize=str(markersize))
        # mrca = tree_phy.common_ancestor({"name": label})
        # mrca.branch_length = mrca.branch_length + 0.1
        # mrca.color = color.split(':')[1]
        nstyle = NodeStyle()
        nstyle["fgcolor"] = color.split(':')[1]
        tree_phy.get_leaves_by_name(label)[0].set_style(nstyle)
    # X = itertools.repeat(float(x), 1)
    # Y = itertools.repeat(float(y), 1)
    # for label, xpt, ypt in zip(labels, X, Y):
    #    plt.text(xpt, ypt, label)

    # for i, (newX, newY) in enumerate(zip(X, Y), start=1):
    #    ax.annotate(city, (newX, newY), xytext=(5, 5), textcoords='offset points', color='g')
    '''for p in place_list:
        place = geo.geocode(city, timeout=2)
        time.sleep(0.5)
        x, y = myMap(p.longitude, p.latitude)
        myMap.plot(x, y, color='g', marker='o', markersize='15')'''
    # Phylo.draw(tree_phy)
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
