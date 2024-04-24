# Version 0.7.0

- Added `tokio` feature for async input/output.

# Version 0.6.0

- Update to rust 2021.
- Better error structure using `thiserror`.
- Use `nom`'s integer parsers

# Version 0.5.1

- Don't panic if `finish` fails while dropping a `Writer`.

# Version 0.5.0

- Re-export the most common structs at crate level.
- Support heavy ion event information.
- Add proper types for energy and length units.

# Version 0.4.1

- Accept (and ignore) empty lines in input.

# Version 0.4.0

- Tweak performance.

# Version 0.3.0

- Try to ensure that the `HepMC` footer is written.

# Version 0.2.1

- Skip lines starting with `HepMC`.

# Version 0.2.0

- `into_inner` methods to retrieve underlying readers and writers.
- Changed behaviour of `Reader::new`.
