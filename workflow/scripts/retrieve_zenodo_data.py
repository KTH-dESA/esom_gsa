import logging, zipfile, requests
import sys

def download_zipped_zenodo_deposit(deposit_url: str, save_path: str):
    try:
        response = requests.get(deposit_url, stream=True)
        if response.status_code == 200:
            with open(save_path, 'wb') as f:
                for chunk in response.iter_content(chunk_size=8192):
                    f.write(chunk)
            print("Download completed.")
        else:
            print(f"Error: Unable to download. Status code: {response.status_code}")
    except requests.exceptions.RequestException as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    logger = logging.getLogger(__name__)

    zenodo_deposit = sys.argv[1]
    save_path = sys.argv[2]
    
    if not save_path[-4:] == ".zip":
        save_path = f"{save_path}.zip"
    
    download_zipped_zenodo_deposit(zenodo_deposit, save_path)

    with zipfile.ZipFile(save_path, "r") as zip_ref:
        zip_ref.extractall(save_path[:-4]) # remove zip extension