import logging, zipfile, requests
from pathlib import Path
from typing import List

def download_zenodo_deposit(doi: str, save_path: str) -> List[str]:
    zenodo_api = "zenodo.org/api/records" 
    deposit_url = Path(zenodo_api, doi)
    
    try:
        response = requests.get(f"https://{str(deposit_url)}")
        if response.status_code == 200:
            files = response.json()["files"] # get file URLs
            
            # download each file
            file_names = []
            for file in files:
                name = Path(save_path, file["key"])
                url = file["links"]["self"]
                
                with requests.get(url, stream=True) as r:
                    with open(name, "wb") as f:
                        for chunk in r.iter_content(chunk_size=8192):
                            if chunk:
                                f.write(chunk)
                print(f"Dowloaded {str(name.absolute())}")
                file_names.append(file["key"])
            print("Download completed")
            return file_names
        else:
            print(f"Error: Unable to download. Status code: {response.status_code}")
            return []
    except requests.exceptions.RequestException as e:
        print(f"Error: {e}")
        return []

if __name__ == "__main__":
    logger = logging.getLogger(__name__)
    
    # zenodo_deposit = "https://zenodo.org/records/10107849"
    zenodo_deposit = snakemake.params.doi
    save_path = snakemake.params.save
    
    files = download_zenodo_deposit(zenodo_deposit, save_path)

    for f in files:
        file_path = Path(save_path, f)
        if str(file_path)[-4:] == ".zip":
            with zipfile.ZipFile(file_path, "r") as zip_ref:
                zip_ref.extractall(str(file_path)[:-4]) # remove zip extension
            file_path.unlink()