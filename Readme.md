# hepmc2

Read and write event files in the `hepmc2` format, also known as
`IO_GenEvent`.

## Caveats

This crate is inspired by the code for the `ReaderAsciiHepMC2` in the
[HepMC3 library](https://gitlab.cern.ch/hepmc/HepMC), version
3.2.0. When using the current version, be aware of

- Lack of rigorous tests
- No support for heavy ions

## Example

```rust
// Read events from `events_in.hepmc2` and write them to `events_out.hepmc2`
use hepmc2::{Reader, Writer};

use std::io::BufReader;
use std::fs::File;

let input = BufReader::new(File::open("events_in.hepmc2")?);
let in_events = Reader::from(input);

let output = File::create("events_out.hepmc2")?;
let mut writer = Writer::try_from(output)?;

for event in in_events {
   let event = event?;
   println!("Current cross section: {}",  event.xs);
   writer.write(&event)?
}
writer.finish()?;
```

License: GPL-3.0-or-later
