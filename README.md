# Markov-Chain Based Procedural Melody Generator
You need Pure Data installed to build or run: [https://puredata.info/downloads/pure-data](https://puredata.info/downloads/pure-data).

To build:
1. Download the [`pd-lib-builder`](https://github.com/pure-data/pd-lib-builder).
2. Place the `pd-lib-builder` in the project repository.
3. Set the `PDINCLUDEDIR` argument in the `Makefile.pdlibbuilder` to `[your pd installation root]/src`

To run:
1. Open `markov_demo.pd` in Pure Data.
2. Turn DSP on.
3. Click a Message block to send melody sequence in (or write one of your own).
4. Click the Bang object indicated as `Start/Stop`. 
