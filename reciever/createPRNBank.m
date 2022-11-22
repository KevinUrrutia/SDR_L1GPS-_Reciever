function [PRN] = createPRNBank()
    PRN = cell([1,34]);
    for i = 1:34
        PRN{i} = generateGoldSeq(i);
    end
end
