import argparse
import os.path

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--file', '-f', type=str, default=None,
                        dest='file', required=True,
                        help='A text file containing region stats')

    args = parser.parse_args()

    # check that file exists

    if not os.path.isfile(args.file):
        print("The file {} does not exist".format(args.file))
        return 1

    # open file, get signal and values
    data = {}
    with open(args.file) as file:
        lines = file.readlines()
        if len(lines) != 9:
            print("Input file has wrong number of lines, please check format")
            return 1

        for lineIndex in range(len(lines)):
            if lineIndex == 0:
                signal = lines[lineIndex].rstrip()
            elif lineIndex == 1:
                data["s_sqrt_spb"] = float(lines[lineIndex])
            elif lineIndex == 2:
                data["sumFullHist"] = float(lines[lineIndex])
            elif lineIndex == 3:
                data["sumFullHistRegion"] = float(lines[lineIndex])
            elif lineIndex == 4:
                data["sumSignalHistRegion"] = float(lines[lineIndex])
            elif lineIndex == 5:
                data["sumFullHistDirect"] = float(lines[lineIndex])
            elif lineIndex == 6:
                data["sumSignalHistDirect"] = float(lines[lineIndex])
            elif lineIndex == 7:
                data["sumFullHistReflected"] = float(lines[lineIndex])
            elif lineIndex == 8:
                data["sumSignalHistReflected"] = float(lines[lineIndex])

    # calculate FOMs

    FOMs = {}
    FOMs['Signal'] = {}
    FOMs['Full'] = {}
    for dataType in ['Signal', 'Full']:
        for numerator in ['sum{}HistRegion'.format(dataType), 'sum{}HistDirect'.format(dataType), 'sum{}HistReflected'.format(dataType)]:
            for denominator in ['sum{}HistRegion'.format(dataType), 'sum{}HistDirect'.format(dataType), 'sum{}HistReflected'.format(dataType)]:
                if denominator != numerator:
                    FOMs[dataType]['{} / {}'.format(numerator, denominator)] = data[numerator] / data[denominator]
                if denominator == 'sum{}HistRegion'.format(dataType) or denominator == 'sum{}HistDirect'.format(dataType):
                    FOMs[dataType]['{} / {} + {}'.format(numerator, denominator, 'sum{}HistReflected'.format(dataType))] = data[numerator] / (data[denominator] + data['sum{}HistReflected'.format(dataType)])
                    if denominator == 'sum{}HistRegion'.format(dataType):
                        FOMs[dataType]['{} / {} + {}'.format(numerator, denominator, 'sum{}HistDirect'.format(dataType))] = data[numerator] / (data[denominator] + data['sum{}HistDirect'.format(dataType)])
                        FOMs[dataType]['{} / {} + {} + {}'.format(numerator, denominator, 'sum{}HistDirect'.format(dataType), 'sum{}HistReflected'.format(dataType))] = data[numerator] / (data[denominator] + data['sum{}HistReflected'.format(dataType)] + data['sum{}HistDirect'.format(dataType)])

    # write FOMs to txt file

    outputFileString = 'FOMs_' + args.file 
    with open(outputFileString, 'w+') as file:
        for dataType in FOMs.keys():
            if dataType == 'Signal':
                file.write('{} light only (MC only):\n'.format(signal))
            else:
                file.write('All PMT hits:\n')
            for FOM in FOMs[dataType].keys():
                file.write('    {}: {}\n'.format(FOM, FOMs[dataType][FOM]))
            file.write('\n')

    return 0

if __name__ == '__main__':
    main()