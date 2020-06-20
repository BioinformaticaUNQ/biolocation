import itertools
import os
import sys

import Bio
from Bio import Entrez
from Bio import SeqIO
from Bio.Alphabet import generic_protein
from geopy.geocoders import Nominatim
import time

home_dir = os.path.expanduser("~")
my_module = os.path.join(home_dir, "Anaconda3\Library\share\lib")
os.environ['PROJ_LIB'] = my_module

from mpl_toolkits.basemap.test import Basemap
import matplotlib.pyplot as plt


def search_features(name, key, table):
    return [element for element in table if element[key] == name]


def search_taxonID(table):
    return [element for element in table if
            element["GBQualifier_name"] == "db_xref" and element["GBQualifier_value"].split(":")[0] == 'taxon']


def fasta_to_genbank(filename):
    file = open(filename, "r")
    # records = [seqrec for seqrec in Bio.SeqIO.parse(file, "fasta", alphabet=generic_protein)]
    records = Bio.SeqIO.parse(file, "fasta", alphabet=generic_protein)
    gb_file = "THIS_IS_YOUR_OUTPUT_FILE.genbank"
    Bio.SeqIO.write(records, gb_file, "genbank")
    file.close()
    return records, gb_file

    # format it for R, e.g.:
    # os.system("sed 's/\t \t/\tNA\t/g; s/\t $/\tNA/g; s/\t|//g' tax_report.txt > all_mammals_WR_NCBI_taxID.txt")
    # taxidTable = pandas.read_table('tax_report.txt')
    # del taxidTable['Numero']
    # taxlist = taxidTable
    # taxlist = pandas.read_table('tax_report.txt', header=0, sep='\t')
    # with open('all_mammals.taxlist.NCBI.txt', "w") as fd:
    #    fd.writelines(str(taxlist))
    # write .table(taxlist, file='all_mammals.taxlist.NCBI.txt', sep='\t', quote=FALSE, row.names = F)
    # taxIdList = [t for t in taxidTable['taxid']]
    # print(taxIdList)
    # return taxIdList


def listIds(seq_record):
    variant = seq_record.description.split('[')[1].replace(']', '')
    Entrez.email = 'myemail@ncbi.nlm.nih.gov'  # provide your email address
    db = 'nucleotide'  # set search to dbVar database
    paramEutils = {'usehistory': 'Y', 'RetMax': '10'}  # Use Entrez search history to cache results
    # download list of GIs for the species
    print('downloading sequences for species organism:' + variant)
    # generate query to Entrez eSearch
    eSearch = Entrez.esearch(db=db, term='("' + variant + '"[Organism])', **paramEutils)
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


def dataset(fileName):
    list_of_records, gb_file = fasta_to_genbank(fileName)
    input_handle = open(gb_file, "r")
    for seq_record in SeqIO.parse(input_handle, "genbank"):
        # for species in list_of_species:
        # trm = 'txid'
        # trm += str(species)
        # trm += '[orgn] AND ddbj embl genbank with limits[filt] NOT transcriptome[All Fields] NOT mRNA[filt] NOT scaffold[All Fields]'
        # build genbank search term..
        # provide email
        # Entrez.email = "paolo_gratton@eva.mpg.de"
        # download list of GIs for the species
        # print('downloading sequences for species txid:' + str(species))
        # handle = Entrez.esearch(db="Nucleotide", term=trm, RetMax="10")
        # record = Entrez.read(handle, validate=False)
        list_of_ids = listIds(seq_record)
        if len(list_of_ids) > 0:
            Entrez.email = "paolo_gratton@eva.mpg.de"  # provide email
            ### read entry in xml and extract features..
            # use Entrez.efetch to read data in .xml format
            handle = Entrez.efetch(db="nucleotide", id=list_of_ids, retmode="xml")
            record = Entrez.read(handle, validate=False)
            print(record)
            # reset lists for values to be stored
            accessions = []
            definitions = []
            organisms = []
            taxonIDs = []
            countries = []
            lat_lons = []
            genes = []
            rRNAs = []
            tRNAs = []
            genera = []
            isolates = []
            haplotypes = []
            pop_variants = []
            dloop = []
            CDSs = []
            organelles = []
            ref_pubmed = []
            ref_authors = []
            ref_journal = []
            ref_titles = []

            feattab = record[0]['GBSeq_feature-table']

            refs = record[0]['GBSeq_references']

            # get and name the 'source' section (s)
            sourcetab = search_features('source', 'GBFeature_key', feattab)
            # output.write(record[i]['GBSeq_primary-accession'] + ": " + str(len(sourcetab))+ "\n") # this was just a check!

            # create temporary objects to collect values of the fields of interest.. note that we run a loop through 'sourcetab'.. It is probably not necessary, since there should be only one 'source' section.. but we do it this way because the 'gene' sections will be often multiple..
            temp_organisms = []
            temp_taxonIDs = []
            temp_countries = []
            temp_lat_lons = []
            temp_isolates = []
            temp_pop_variants = []
            temp_haplotypes = []
            temp_organelles = []

            if len(sourcetab) > 0:
                for j in range(len(sourcetab)):
                    if len(search_features('organism', 'GBQualifier_name', sourcetab[j][
                        'GBFeature_quals'])) > 0:  # store the organism name (Linnean binomial)
                        temp_organisms.append(
                            search_features('organism', 'GBQualifier_name', sourcetab[j]['GBFeature_quals'])[0][
                                'GBQualifier_value'])  # we assume that there is only one field 'organism' per each source field..
                    else:
                        temp_organisms.append('NA')

                    # store organelle information
                    if len(search_features('organelle', 'GBQualifier_name',
                                           sourcetab[j]['GBFeature_quals'])) > 0:
                        temp_organelles.append(
                            search_features('organelle', 'GBQualifier_name', sourcetab[j]['GBFeature_quals'])[
                                0][
                                'GBQualifier_value'])  # we assume that there is only one field 'organelle' per each source field..
                    else:
                        temp_organelles.append('NA')

                        # store the isolate name
                    if len(search_features('isolate', 'GBQualifier_name', sourcetab[j]['GBFeature_quals'])) > 0:
                        temp_isolates.append(
                            search_features('isolate', 'GBQualifier_name', sourcetab[j]['GBFeature_quals'])[0][
                                'GBQualifier_value'])  # we assume that there is only one field 'isolate' per each source field..
                    else:
                        temp_isolates.append('NA')

                    # store the pop_variant
                    if len(search_features('pop_variant', 'GBQualifier_name',
                                           sourcetab[j]['GBFeature_quals'])) > 0:
                        temp_pop_variants.append(
                            search_features('pop_variant', 'GBQualifier_name', sourcetab[j]['GBFeature_quals'])[
                                0][
                                'GBQualifier_value'])  # we assume that there is only one field 'isolate' per each source field..
                    else:
                        temp_pop_variants.append('NA')

                    # store the haplotype
                    if len(search_features('haplotype', 'GBQualifier_name',
                                           sourcetab[j]['GBFeature_quals'])) > 0:
                        temp_haplotypes.append(
                            search_features('haplotype', 'GBQualifier_name', sourcetab[j]['GBFeature_quals'])[
                                0][
                                'GBQualifier_value'])  # we assume that there is only one field 'isolate' per each source field..
                    else:
                        temp_haplotypes.append('NA')

                        # store the NCBI taxonID
                    if len(search_taxonID(sourcetab[0]['GBFeature_quals'])) > 0:
                        temp_taxonIDs.append(
                            search_taxonID(sourcetab[0]['GBFeature_quals'])[0]['GBQualifier_value'].split(':')[
                                1])  # we assume that there is only one sub-field 'taxon:' per each source field..
                    else:
                        temp_taxonIDs.append('NA')

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

                organisms.append(','.join(temp_organisms))
                taxonIDs.append(','.join(temp_taxonIDs))
                countries.append(','.join(temp_countries))
                lat_lons.append(','.join(temp_lat_lons))
                isolates.append(','.join(temp_isolates))
                haplotypes.append(','.join(temp_haplotypes))
                pop_variants.append(','.join(temp_pop_variants))
                organelles.append(','.join(temp_organelles))
            else:
                organisms.append('NA')
                taxonIDs.append('NA')
                countries.append('NA')
                lat_lons.append('NA')
                isolates.append('NA')
                haplotypes.append('NA')
                pop_variants.append('NA')
                organelles.append('NA')
            print('lat_lons')
            print(lat_lons[0])
        draw(countries[0].split(',')[0].split(': ')[1])

        # loop through all records
        '''for i in range(len(record)):
            if float(record[i]['GBSeq_length']) < 20001:  # I am ignoring all sequences longer than 20000 bp!! (they tend to be fucking genomic scaffolds..), but it will likely not be enough to let it go through everything!!!
                # genera.append(genus)  # store the genus for this sequence
                accessions.append(record[i]['GBSeq_primary-accession'])  # store the sequence accession number
                definitions.append(record[i]['GBSeq_definition'])  # store sequence definition (description)
                # name the 'feature table' section in the entry
                feattab = record[i]['GBSeq_feature-table']

                refs = record[i]['GBSeq_references']  # I decided to take only the first reference, at the moment..
                if len(re.findall("GBReference_pubmed", str(record[i]['GBSeq_references'][0].keys()))) > 0:
                    ref_pubmed.append(record[i]['GBSeq_references'][0]['GBReference_pubmed'])
                else:
                    ref_pubmed.append('NA')
                if len(re.findall("GBReference_authors", str(record[i]['GBSeq_references'][0].keys()))) > 0:
                    ref_authors.append(record[i]['GBSeq_references'][0]['GBReference_authors'])
                else:
                    ref_authors.append('NA')
                if len(re.findall("GBReference_title", str(record[i]['GBSeq_references'][0].keys()))) > 0:
                    ref_titles.append(record[i]['GBSeq_references'][0]['GBReference_title'])
                else:
                    ref_titles.append('NA')
                if len(re.findall("GBReference_journal", str(record[i]['GBSeq_references'][0].keys()))) > 0:
                    ref_journal.append(record[i]['GBSeq_references'][0]['GBReference_journal'])
                else:
                    ref_journal.append('NA')

                # get and name the 'source' section (s)
                sourcetab = search_features('source', 'GBFeature_key', feattab)
                # output.write(record[i]['GBSeq_primary-accession'] + ": " + str(len(sourcetab))+ "\n") # this was just a check!

                # create temporary objects to collect values of the fields of interest.. note that we run a loop through 'sourcetab'.. It is probably not necessary, since there should be only one 'source' section.. but we do it this way because the 'gene' sections will be often multiple..
                temp_organisms = []
                temp_taxonIDs = []
                temp_countries = []
                temp_lat_lons = []
                temp_isolates = []
                temp_pop_variants = []
                temp_haplotypes = []
                temp_organelles = []

                if len(sourcetab) > 0:
                    for j in range(len(sourcetab)):
                        if len(search_features('organism', 'GBQualifier_name', sourcetab[j][
                            'GBFeature_quals'])) > 0:  # store the organism name (Linnean binomial)
                            temp_organisms.append(
                                search_features('organism', 'GBQualifier_name', sourcetab[j]['GBFeature_quals'])[0][
                                    'GBQualifier_value'])  # we assume that there is only one field 'organism' per each source field..
                        else:
                            temp_organisms.append('NA')

                        # store organelle information
                        if len(search_features('organelle', 'GBQualifier_name', sourcetab[j]['GBFeature_quals'])) > 0:
                            temp_organelles.append(
                                search_features('organelle', 'GBQualifier_name', sourcetab[j]['GBFeature_quals'])[0][
                                    'GBQualifier_value'])  # we assume that there is only one field 'organelle' per each source field..
                        else:
                            temp_organelles.append('NA')

                            # store the isolate name
                        if len(search_features('isolate', 'GBQualifier_name', sourcetab[j]['GBFeature_quals'])) > 0:
                            temp_isolates.append(
                                search_features('isolate', 'GBQualifier_name', sourcetab[j]['GBFeature_quals'])[0][
                                    'GBQualifier_value'])  # we assume that there is only one field 'isolate' per each source field..
                        else:
                            temp_isolates.append('NA')

                        # store the pop_variant
                        if len(search_features('pop_variant', 'GBQualifier_name', sourcetab[j]['GBFeature_quals'])) > 0:
                            temp_pop_variants.append(
                                search_features('pop_variant', 'GBQualifier_name', sourcetab[j]['GBFeature_quals'])[0][
                                    'GBQualifier_value'])  # we assume that there is only one field 'isolate' per each source field..
                        else:
                            temp_pop_variants.append('NA')

                        # store the haplotype
                        if len(search_features('haplotype', 'GBQualifier_name', sourcetab[j]['GBFeature_quals'])) > 0:
                            temp_haplotypes.append(
                                search_features('haplotype', 'GBQualifier_name', sourcetab[j]['GBFeature_quals'])[0][
                                    'GBQualifier_value'])  # we assume that there is only one field 'isolate' per each source field..
                        else:
                            temp_haplotypes.append('NA')

                            # store the NCBI taxonID
                        if len(search_taxonID(sourcetab[0]['GBFeature_quals'])) > 0:
                            temp_taxonIDs.append(
                                search_taxonID(sourcetab[0]['GBFeature_quals'])[0]['GBQualifier_value'].split(':')[
                                    1])  # we assume that there is only one sub-field 'taxon:' per each source field..
                        else:
                            temp_taxonIDs.append('NA')

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

                    organisms.append(','.join(temp_organisms))
                    taxonIDs.append(','.join(temp_taxonIDs))
                    countries.append(','.join(temp_countries))
                    lat_lons.append(','.join(temp_lat_lons))
                    isolates.append(','.join(temp_isolates))
                    haplotypes.append(','.join(temp_haplotypes))
                    pop_variants.append(','.join(temp_pop_variants))
                    organelles.append(','.join(temp_organelles))
                else:
                    organisms.append('NA')
                    taxonIDs.append('NA')
                    countries.append('NA')
                    lat_lons.append('NA')
                    isolates.append('NA')
                    haplotypes.append('NA')
                    pop_variants.append('NA')
                    organelles.append('NA')
                print('countries')
                print(countries)
                print('lat_lons')
                print(lat_lons)'''


def draw(city):
    place_list = []
    geo = Nominatim(user_agent='BioLocation', timeout=2)
    plt.figure(figsize=(16, 12))
    # fig, ax = plt.subplots()
    myMap = Basemap(projection='robin', lon_0=0, lat_0=0) # ax=ax
    place = geo.geocode(city)
    print(place)
    x, y = myMap(place.longitude, place.latitude)
    myMap.drawcoastlines()
    myMap.drawcountries()
    myMap.fillcontinents(color='brown')
    myMap.drawmapboundary(fill_color='aqua')
    # myMap.bluemarble()
    # ax.text(x, y, 'label', ha='center', size=20)
    # myMap.plot(x, y, color='g', marker='o', markersize='15')
    # ax.scatter(x, y)
    X = itertools.repeat(float(x), 1)
    Y = itertools.repeat(float(y), 1)
    labels = [city]
    for label, xpt, ypt in zip(labels, X, Y):
        plt.text(xpt, ypt, label)
    # for i, (newX, newY) in enumerate(zip(X, Y), start=1):
    #    ax.annotate(city, (newX, newY), xytext=(5, 5), textcoords='offset points', color='g')
    '''for p in place_list:
        place = geo.geocode(city, timeout=2)
        time.sleep(0.5)
        x, y = myMap(p.longitude, p.latitude)
        myMap.plot(x, y, color='g', marker='o', markersize='15')'''
