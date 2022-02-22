"""
Rename pdfs in a folder to the contents title, if the document has a title.
Needs pdfrw to work.
If that module is not installed one is prompted during runtime if one wants to install it.
Can be run directly or by using the fct: rnPdfsByTitle.
If collisions occur, idk what happens, depends on os.rename.
"""
# pdf rename
# requires pdfrw
import os

# import pdfrw and prompt if not installed
# could theoretically be removed as python would tell that the package could not be found
try:
    from pdfrw import PdfReader
except ImportError:
    ans = input("package pdfrw missing. Want to install it? (y/n):")
    if "y" in ans.lower():
        osres = os.popen("pip install pdfrw")
        print("result from trying to install pdfrw:", osres.read())
        osres.close()
        from pdfrw import PdfReader

from pathlib import Path
import re
from tkinter import filedialog


def rnPdfsByTitle(pat, usePrint=True):
    resl = []  # (original filename, title)

    # read all pdf in the current folder and extract the title
    # could prbly be speed up by using multiprocess but not that important.
    for f in pat.iterdir():
        op = pat/f
        if op.suffix.lower() != ".pdf":
            continue
        reader = PdfReader(op)
        # [1:-1] : remove parenthesis
        title = reader.Info.Title
        if title is not None and title != "":
            resl.append((op, reader.Info.Title[1:-1]))
        elif usePrint:
            print("(W): title is none, keep name", op, title)

    # rename the files, make the new filename work with windows os.
    for tpl in resl:
        cleanName = tpl[1].replace(" ", "_")+".pdf"
        cleanName = re.sub(r'[^\w\-_\. ]', '_', cleanName)
        newName = tpl[0].parent/cleanName
        os.rename(str(tpl[0]), str(newName))

if __name__ == "__main__":
    print("This script tries to rename pdf files based on title (inplace)")
    p = filedialog.askdirectory()
    pat = Path(p)
    print("Reanaming files with titles in: ", p)
    rnPdfsByTitle(pat)
    # if run directly, wait for input to close the console
    input()
