import numpy as np

def get_Binary(data):
    binary = []
    for d in data:
        # I_data.append(decimalToBinary(data))
        binary.append(format(d,'08b'))

    # print('Binary', binary)

    return binary


def decimalToBinary(n):
    return bin(n).replace("0b", "")

def save_to_binary(data, file):

    f = open(file, 'w+b')
    byte_arr = data
    binary_format = bytearray(byte_arr)
    f.write(binary_format)
    f.close()

def get_number(sz):
    data = []

    for i in range(sz):
        temp = np.random.choice(10,1) 
        # print(type(temp[0]))
        a = temp[0].astype(np.uint8)
        data.append(a)
    
    # print('Decimal: ', data)

    binary = get_Binary(data)

    # print('Binary: ', binary)

    return data

def get_IQ_Data(I,Q):

    """
    Creates IQ Data and saves to binary IQ Data File

    Args:
        I: Psuedo-I data for IQ Data
        Q: Pseudo-Q data for IQ Data

    """

    IQ = []

    for i in range(len(I)): # Does alternate append for I-Q-I-Q Sequence
        IQ.append(I[i])
        IQ.append(Q[i])

    # print('Decimal: ', IQ)

    IQ_Binary = get_Binary(IQ)
    # print('Binary: ', IQ_Binary)

    save_to_binary(IQ, 'IQ_Data') # Saves IQ_Data file 



if __name__ == '__main__':

    size = 1023

    I = get_number(size)
    Q = get_number(size)

    get_IQ_Data(I,Q)
    
    print('Psuedo IQ Data Saved')
    

