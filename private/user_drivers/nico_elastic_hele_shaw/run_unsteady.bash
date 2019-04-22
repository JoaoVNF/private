# full simulation
./elastic_hele_shaw_displ_unsteady --dt 0.01 --epsilon_t 0.1 --l_bubble 0.0 --pre_stress 30.0e3 --r_bubble 0.3 --adapt_frequency 1 --dir RESLT/ &

# Foeppl-von-Karman for the elastic sheet only
./elastic_hele_shaw_displ_unsteady --pin_hs --pre_stress 30.0e3 --adapt_frequency 1 --dir RESLT/ &