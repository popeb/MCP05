

ANCHOR_RAD = 500000
PEAKS      = $(wildcard peaks/*.bed)
BIN        = 200
ARCHNAME   = ../W$(ANCHOR_RAD)_B$(BIN)
TRSIZE     = 200
NTRAINS    = 10

.PHONY : train

train : train_dae.cfg $(ARCHNAME).h5 
	@echo "training ${NTRAINS} times ... it's gonna take a while (a couple of days on my machine)"
	for i in {1..${NTRAINS}}; do dimer_train_dae.py --raw $(ARCHNAME).h5:train$(TRSIZE) $< $(shell ./hiddens.py --nin 7 --nout 5 $(ARCHNAME).h5 4); done

