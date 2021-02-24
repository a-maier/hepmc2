# hepmc2

Read and write event files in the `hepmc2` format, also known as
`IO_GenEvent`.

# Caveats

This crate is inspired by the code for the `ReaderAsciiHepMC2` in the
[HepMC3 library](https://gitlab.cern.ch/hepmc/HepMC), version
3.2.0. The current version is not ready for production use. In
particular, be aware of

- Non-existing documentation
- Lack of rigorous tests
- No support for heavy ions

# Examples

```rust,no_run
// Read events from `events_in.hepmc2` and write them to `events_out.hepmc2`
use hepmc2::reader::Reader;
use hepmc2::writer::Writer;

use std::io::BufReader;
use std::fs::File;

let input = BufReader::new(File::open("events_in.hepmc2")?);
let in_events = Reader::from(input);

let output = File::create("events_out.hepmc2")?;
let mut writer = Writer::try_from(output)?;

for event in in_events {
   let event = event?;
   writer.write(&event)?
}
writer.finish()?
```
