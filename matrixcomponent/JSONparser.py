import json

# class Parser(object):
#
#     infile = ''
#
#     def __init__(self, file):
#         self.infile = file

    def parse(self, file):
        with open(file, "r") as read_file:
            data = json.load(read_file)


