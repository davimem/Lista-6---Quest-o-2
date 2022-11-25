set terminal pdf
set output "Diferenca.pdf" 
plot "saida.dat" w lp pt 7 t "Abs(Numerico-Exato)"
set terminal x11
rep