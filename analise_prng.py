import numpy as np
import sys
import math

numeros = []
maior = 0
menor = sys.maxsize * 2 + 1
with open("numeros_mt.txt") as file:
    num_str = file.readline()
    while num_str:
        num = int(num_str)
        if num > maior:
            maior = num
        if num < menor:
            menor = num
        numeros.append(num)
        num_str = file.readline()

media_esp = (maior + menor) / 2
dp_esp = math.pow(( math.pow((maior-menor),2) / 12 ),(1/2))
media_real = np.average(numeros)
dp_real = np.std(numeros)

print("Maior número: ", maior)
print("Menor número: ", menor)
print("Média esperada: ", media_esp)
print("Média real: ", media_real)
print("Desvio padrão esperado: ", dp_esp)
print("Desvio padrão real: ", dp_real)