% function binaryVector = convertToBinaryVector(inputVector)
% 
%     binaryVector = []; 
% 
%     for i = 1:length(inputVector)
%         binArray = dec2bin(inputVector(i), 2) - '0'; 
%         binaryVector = [binaryVector, binArray];
%     end
% end
function binaryVector = convertToBinaryVector(inputVector, M)
    bitLength = ceil(log2(M));
    binaryVector = []; 
    for i = 1:length(inputVector)
        binArray = dec2bin(inputVector(i), bitLength) - '0'; 
        binaryVector = [binaryVector, binArray];
    end
end
