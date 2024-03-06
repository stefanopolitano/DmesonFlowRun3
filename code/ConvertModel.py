'''
python script to convert ML models from .pkl to .model
'''

import os
import argparse

from hipe4ml.model_handler import ModelHandler
from hipe4ml_converter.h4ml_converter import H4MLConverter

parser = argparse.ArgumentParser(description='Arguments to pass')
parser.add_argument('inFilePkl', metavar='text', default='model.pkl', help='input pickle file to be converted')
parser.add_argument('nfeatures', metavar='int', default=10, help='number of features')
args = parser.parse_args()

ModelPath = os.path.expanduser(args.inFilePkl)
print(f'Loaded saved model: {ModelPath}')
ModelHandl = ModelHandler()
ModelHandl.load_model_handler(ModelPath)
model_converter = H4MLConverter(ModelHandl) # create the converter object
model_onnx = model_converter.convert_model_onnx(1, args.nfeatures) # convert the model to ONNX format

outFileName = ModelPath.replace('.pkl', '.onnx')

model_converter.dump_model_onnx(model_onnx, outFileName) # save the model to a file
print(f'Saved model: {outFileName}')