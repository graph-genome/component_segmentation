import os
from glob import glob

if __name__ == "__main__":

    files = glob(r'D:\josiah\Projects\component_segmentation\matrixcomponent\data\lung_cancer.w1000\*schematic.json')
    for f in files:
        print(f)
        contents = open(f, 'r').read()
        with open('data/shortened/' + os.path.basename(f), 'w') as out:
            for c in contents:
                if c != ' ' and c != '\n':
                    out.write(c)
