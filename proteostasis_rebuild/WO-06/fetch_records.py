"""WO-06: fetch the bibliographic records every parameter citation is checked against.

two kinds of record are saved under records/, never edited by hand:

  pubmed_<pmid>.xml   NCBI efetch (db=pubmed) record: journal, year, volume,
                      pages, authors, abstract. this is the bibliographic check.
  pmc_<pmcid>.xml     NCBI efetch (db=pmc) full text, when an open-access PMC
                      copy exists. used only to locate a quoted value.
  citmatch.json       NCBI ecitmatch result for each citation AS THE LEGACY
                      WROTE IT (journal|year|volume|first page). NOT_FOUND here
                      is the machine evidence behind a MISCITED/UNVERIFIED label.

run once with network access; the audit and tests then read only records/.
"""
import json
import re
import time
import urllib.parse
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
REC = HERE / "records"
E = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"

# every pmid the audit relies on. keep in sync with parameter_audit.tsv
# (test_wo06 checks that each pmid cited in the tsv has a saved record)
PMIDS = [
    "16916930",  # Belle 2006 PNAS, yeast half-lives
    "25466257",  # Christiano 2014 Cell Rep, yeast turnover
    "9223639",   # Pierpaoli 1997 J Mol Biol, DnaK power stroke
    "9843444",   # Pierpaoli 1998 Biochemistry, DnaK/DnaJ peptide rates
    "9506960",   # Pierpaoli 1998 JBC, substoichiometric DnaJ/GrpE
    "8566548",   # Lorimer 1996 FASEB J, chaperonin quantitative assessment
    "22479486",  # Upadhyay 2012 PLoS One, inclusion-body kinetics
    "11135201",  # Hoffmann, Posten, Rinas 2001 Biotechnol Bioeng
    "18662548",  # Drummond & Wilke 2008 Cell
    "19763154",  # Drummond & Wilke 2009 Nat Rev Genet
    "10601016",  # Mogk 1999 EMBO J, thermolabile proteins, DnaK/ClpB
    "24183671",  # Ciryam 2013 Cell Rep, supersaturation
    "23256155",  # Ciryam 2013 PNAS 110:E132 (the only 2013 PNAS by Ciryam)
    "23894132",  # Bednarska 2013 Microbiology
    "24239291",  # what Mol Cell 52:617 (2013) actually is
    "24867638",  # De Los Rios & Barducci 2014 eLife
    "38421032",  # Landerer 2024 MBE
    "42406629",  # Stikeleather 2026 NAR
    "24114984",  # Milo 2013 BioEssays, protein molecules per cell volume
    "22832197",  # Calloni 2012 Cell Rep, DnaK hub
    "26641532",  # Schmidt 2016 Nat Biotechnol, condition-dependent proteome
    "24766808",  # Li 2014 Cell, absolute synthesis rates
    "4551144",   # Goldberg 1972 PNAS, degradation of abnormal proteins in E. coli
    "6989832",   # Larrabee 1980 JBC, synthesis vs degradation in growing E. coli
    "4912536",   # Nath & Koch 1970 JBC, rapidly/slowly decaying protein components
    "21829590",  # Volkmer & Heinemann 2011 PLoS One, cell volume (Schmidt's volumes)
]

# citations exactly as the legacy wrote them: journal|year|volume|first_page
CITMATCH = {
    "Ciryam2013_PNAS_110_E3453": "proc natl acad sci u s a|2013|110|E3453",
    "Bednarska2013_MolCell_52_617": "mol cell|2013|52|617",
    "Stirling2018_CellRep_25_2242": "cell rep|2018|25|2242",
    "DrummondWilke2009_Cell": "cell|2009||",
    "Yamanaka2017_CurrBiol": "curr biol|2017||",
    "Pierpaoli1997_EMBOJ": "embo j|1997||",
    "Mogk1999_EMBOJ_18": "embo j|1999|18|6934",
    "Lorimer1996_FASEBJ_10": "faseb j|1996|10|5",
    "Upadhyay2012_PLoSOne_7": "plos one|2012|7|e33951",
    "Belle2006_PNAS_103": "proc natl acad sci u s a|2006|103|13004",
}


def get(url, tries=4):
    for i in range(tries):
        try:
            return urllib.request.urlopen(url, timeout=60).read().decode()
        except Exception:
            time.sleep(2 * (i + 1))
    raise RuntimeError(url)


def efetch(db, uid):
    q = urllib.parse.urlencode({"db": db, "id": uid, "retmode": "xml"})
    time.sleep(0.4)
    return get(E + "efetch.fcgi?" + q)


def main():
    REC.mkdir(exist_ok=True)
    for pmid in PMIDS:
        f = REC / f"pubmed_{pmid}.xml"
        if not f.exists():
            f.write_text(efetch("pubmed", pmid))
        x = f.read_text()
        m = re.search(r'<ArticleId IdType="pmc">(PMC\d+)</ArticleId>', x.split("<ReferenceList")[0])
        if m:
            g = REC / f"pmc_{m.group(1)}.xml"
            if not g.exists():
                g.write_text(efetch("pmc", m.group(1)[3:]))
    cm = {}
    for key, bdata in CITMATCH.items():
        q = urllib.parse.urlencode({"db": "pubmed", "retmode": "xml", "bdata": bdata + f"|x|{key}|"})
        time.sleep(0.4)
        cm[key] = {"query": bdata, "result": get(E + "ecitmatch.cgi?" + q).strip().split("|")[-1]}
    (REC / "citmatch.json").write_text(json.dumps(cm, indent=2))
    print(json.dumps(cm, indent=1))


if __name__ == "__main__":
    main()
