import selenium
from selenium import webdriver
from selenium.webdriver.common.keys import Keys
import re
import requests
import time
import pandas as pd
results = []

def start_job(pdb_id,chain, results):
    driver = webdriver.Firefox()
    driver.get("http://consurf.tau.ac.il/index_from_ConSurfDB.php?chain="+chain+"&pdb_ID="+pdb_id)
    time.sleep(10)
    elem = driver.find_element_by_name("submitForm")
    elem.click()
    time.sleep(10)
    job = re.findall("ConSurf run no. (\d*) PDB ID: (.*), Chain: (\w) <" , driver.page_source)
    job_id = None
    try:
        job_id = job[0][0]
        pdb_id = job[0][1]
        chain = job[0][2]
    except Exception as e:
        print(e)
    results.append({"job_id":job_id,"pdb_id":pdb_id,"chain":chain})
    driver.close()
    return results

#for i in codes:
#    results = []
#    pdb_id = codes[0]
#    chain = codes[1]
#    results = start_job(pdb_id,chain,results)
#    df = pd.DataFrame(results)
#    df.to_csv(path, index=False)


# response = requests.get("http://consurf.tau.ac.il/results/job_id/output.php".replace("job_id", job_id))
