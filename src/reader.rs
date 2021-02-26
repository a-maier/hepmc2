use std::collections::BTreeMap;
use std::default::Default;
use std::fmt::{self, Display};
use std::io::{self, BufRead};
use std::iter::Iterator;
use std::num::{ParseFloatError, ParseIntError};

use crate::event::*;

use nom::{
    bytes::complete::{take_until, take_while1},
    character::complete::{char, space1},
    combinator::{opt, recognize},
    number::complete::{double, recognize_float},
    sequence::{delimited, preceded, tuple},
    IResult,
};

const BUF_SIZE: usize = 256;

/// Reader for the HepMC2 format
#[derive(Debug, Eq, PartialEq, Ord, PartialOrd, Hash, Default)]
pub struct Reader<T> {
    stream: T,
    line: String,
    line_nr: usize,
}

impl<T> Reader<T> {
    /// Retrieve the underlying reader
    pub fn into_inner(self) -> T {
        self.stream
    }
}

impl<T: BufRead> Reader<T> {
    pub fn new(stream: T) -> Self {
        stream.into()
    }
}

impl<T: BufRead> From<T> for Reader<T> {
    fn from(stream: T) -> Self {
        Self {
            stream,
            line: String::with_capacity(BUF_SIZE),
            line_nr: 0,
        }
    }
}

impl<T: BufRead> Reader<T> {
    fn skip_headers(&mut self) -> Result<(), io::Error> {
        while self.line.trim().is_empty() || self.line.starts_with("HepMC") {
            self.line.clear();
            if self.stream.read_line(&mut self.line)? == 0 {
                break;
            }
            self.line_nr += 1;
        }
        Ok(())
    }

    fn parse_event_inner(&mut self) -> Result<Event, ParseError> {
        let mut event = parse_event_line(&self.line)?;
        loop {
            self.line.clear();
            if self.stream.read_line(&mut self.line)? == 0 {
                break;
            };
            self.line_nr += 1;
            match self.line.as_bytes().first() {
                Some(b'E') => break,
                Some(b'V') => parse_vertex_line(&self.line, &mut event)?,
                Some(b'P') => parse_particle_line(&self.line, &mut event)?,
                Some(b'U') => parse_units_line(&self.line, &mut event)?,
                Some(b'F') => parse_pdf_info_line(&self.line, &mut event)?,
                Some(b'H') => {
                    if self.line.starts_with("HepMC") { continue; }
                    parse_heavy_ion_line(&self.line, &mut event)?
                },
                Some(b'N') => parse_weight_names_line(&self.line, &mut event)?,
                Some(b'C') => parse_xs_info_line(&self.line, &mut event)?,
                _ => return Err(ParseError::BadPrefix),
            };
        }
        Ok(event)
    }

    fn parse_event(&mut self) -> Result<Event, LineParseError> {
        self.parse_event_inner().map_err(|err| LineParseError {
            err,
            line: self.line.clone(),
            line_nr: self.line_nr,
        })
    }
}

fn whitespace(line: &str) -> IResult<&str, &str> {
    space1(line)
}

fn non_whitespace(line: &str) -> IResult<&str, &str> {
    take_while1(|c: char| !c.is_ascii_whitespace())(line)
}

fn int(line: &str) -> IResult<&str, &str> {
    let mut int_parser = recognize(tuple((
        opt(char('-')),
        take_while1(|c: char| c.is_ascii_digit()),
    )));
    int_parser(line)
}

fn ws_int(line: &str) -> IResult<&str, &str> {
    preceded(whitespace, int)(line)
}

fn ws_double(line: &str) -> IResult<&str, f64> {
    preceded(whitespace, double)(line)
}

fn string(line: &str) -> IResult<&str, &str> {
    delimited(char('"'), take_until("\""), char('"'))(line)
}

fn parse_event_line(line: &str) -> Result<Event, ParseError> {
    let (rest, parsed) = tuple((
        char('E'),
        ws_int, // event number
        ws_int, // mpi
        whitespace,
        recognize_float,
        whitespace, //event scale
        recognize_float,// alpha_qcd
        whitespace,
        recognize_float,// alpha_qed
        ws_int, // signal_process_id
        ws_int, // signal_process_vertex
        ws_int, // num_vertices
        ws_int, // beam1
        // have to stop here, tuple can't handle more entries
    ))(line)?;
    let (
        _,
        event_number,
        mpi,
        _,
        event_scale,
        _,
        alpha_qcd,
        _,
        alpha_qed,
        signal_process_id,
        signal_process_vertex,
        num_vertices,
        _beam1,
    ) = parsed;
    let (mut rest, parsed) = tuple((
        ws_int, // beam2
        ws_int, // random states
    ))(rest)?;
    let (_beam2, nrandom_states) = parsed;
    let nrandom_states = nrandom_states.parse()?;
    let mut random_states = Vec::with_capacity(nrandom_states);
    for _ in 0..nrandom_states {
        let (rem, random_state) = preceded(whitespace, int)(rest)?;
        rest = rem;
        let random_state = random_state.parse()?;
        random_states.push(random_state);
    }
    let (mut rest, parsed) = preceded(whitespace, int)(rest)?;
    let nweights = parsed.parse()?;
    let mut weights = Vec::with_capacity(nweights);
    for _ in 0..nweights {
        let (rem, weight) = ws_double(rest)?;
        rest = rem;
        weights.push(weight);
    }
    let event = Event {
        number: event_number.parse()?,
        mpi: mpi.parse()?,
        scale: event_scale.parse()?,
        alpha_qcd: alpha_qcd.parse()?,
        alpha_qed: alpha_qed.parse()?,
        signal_process_id: signal_process_id.parse()?,
        signal_process_vertex: signal_process_vertex.parse()?,
        random_states,
        weights,
        vertices: Vec::with_capacity(num_vertices.parse()?),
        weight_names: Default::default(),
        xs: Default::default(),
        energy_unit: Default::default(),
        length_unit: Default::default(),
        pdf_info: Default::default(),
    };
    Ok(event)
}

fn parse_vertex_line(line: &str, event: &mut Event) -> Result<(), ParseError> {
    let (mut rest, parsed) = tuple((
        char('V'),
        ws_int,     // barcode
        ws_int,     // status
        whitespace,
        recognize_float, // x
        whitespace,
        recognize_float, // y
        whitespace,
        recognize_float, // z
        whitespace,
        recognize_float, // t
        ws_int,        // num_orphans_in
        ws_int,        // num_particles_out
        ws_int,        // num_weights
    ))(line)?;
    let (
        _,
        barcode,
        status,
        _,
        x,
        _,
        y,
        _,
        z,
        _,
        t,
        _num_orphans_in,
        num_particles_out,
        num_weights,
    ) = parsed;
    let num_weights = num_weights.parse()?;
    let mut weights = Vec::with_capacity(num_weights);
    for _ in 0..num_weights {
        let (rem, weight) = ws_double(rest)?;
        rest = rem;
        weights.push(weight);
    }
    let vertex = Vertex {
        barcode: barcode.parse()?,
        status: status.parse()?,
        x: x.parse()?,
        y: y.parse()?,
        z: z.parse()?,
        t: t.parse()?,
        weights,
        particles_in: Vec::new(),
        particles_out: Vec::with_capacity(num_particles_out.parse()?),
    };
    event.vertices.push(vertex);
    Ok(())
}

fn parse_particle_line(
    line: &str,
    event: &mut Event,
) -> Result<(), ParseError> {
    let (rest, parsed) = tuple((
        char('P'),
        ws_int, // barcode
        ws_int, // id
        whitespace,
        recognize_float,// px
        whitespace,
        recognize_float,// py
        whitespace,
        recognize_float,// pz
        whitespace,
        recognize_float,// E
        whitespace,
        recognize_float,// m
        whitespace,
    ))(line)?;
    let (_, _barcode, id, _, px, _, py, _, pz, _, e, _, m, _) = parsed;
    let (mut rest, parsed) = tuple((
        int,// status
        whitespace,
        recognize_float,// theta
        whitespace,
        recognize_float,// phi
        ws_int, // end_vtx_code
        ws_int        // flowsize
    ))(rest)?;
    let (status, _, theta, _, phi, end_vtx_code, flowsize) = parsed;
    let mut flows = BTreeMap::new();
    for _ in 0..flowsize.parse()? {
        let (rem, flowidx) = preceded(whitespace, int)(rest)?;
        let flowidx = flowidx.parse()?;
        let (rem, flowval) = preceded(whitespace, int)(rem)?;
        let flowval = flowval.parse()?;
        rest = rem;
        flows.insert(flowidx, flowval);
    }
    let particle = Particle {
        id: id.parse()?,
        p: FourVector::txyz(e.parse()?, px.parse()?, py.parse()?, pz.parse()?),
        m: m.parse()?,
        status: status.parse()?,
        theta: theta.parse()?,
        phi: phi.parse()?,
        flows,
        end_vtx: end_vtx_code.parse()?,
    };
    // TODO: handling of end_vtx is ReaderAsciiHepMC2.cc is obscure and undocumented
    if let Some(vertex) = event.vertices.last_mut() {
        if particle.end_vtx == vertex.barcode {
            vertex.particles_in.push(particle);
        } else {
            vertex.particles_out.push(particle);
        }
    } else {
        return Err(ParseError::NoVertex);
    }
    Ok(())
}

fn parse_units_line(line: &str, event: &mut Event) -> Result<(), ParseError> {
    let (_rest, parsed) = tuple((
        char('U'),
        whitespace,
        non_whitespace,
        whitespace,
        non_whitespace,
    ))(line)?;
    let (_, _, energy, _, length) = parsed;
    event.energy_unit = energy.to_owned();
    event.length_unit = length.to_owned();
    Ok(())
}

fn parse_pdf_info_line(
    line: &str,
    event: &mut Event,
) -> Result<(), ParseError> {
    let (_rest, parsed) = tuple((
        char('F'),
        ws_int, // id0
        ws_int, // id1
        whitespace,
        recognize_float,// x0
        whitespace,
        recognize_float,// x1
        whitespace,
        recognize_float,// scale
        whitespace,
        recognize_float,// xf0
        whitespace,
        recognize_float,// xf1
        whitespace,
        opt(int),// pdf_id0
        whitespace,
        opt(int),   // pdf_id1
    ))(line)?;
    let (
        _,
        id0,
        id1,
        _,
        x0,
        _,
        x1,
        _,
        scale,
        _,
        xf0,
        _,
        xf1,
        _,
        pdf_id0,
        _,
        pdf_id1,
    ) = parsed;
    let pdf_info = PdfInfo {
        parton_id: [id0.parse()?, id1.parse()?],
        x: [x0.parse()?, x1.parse()?],
        scale: scale.parse()?,
        xf: [xf0.parse()?, xf1.parse()?],
        pdf_id: [
            pdf_id0.map(|id| id.parse()).transpose()?.unwrap_or(0),
            pdf_id1.map(|id| id.parse()).transpose()?.unwrap_or(0),
        ],
    };
    event.pdf_info = pdf_info;
    Ok(())
}

fn parse_heavy_ion_line(
    _line: &str,
    _event: &mut Event,
) -> Result<(), ParseError> {
    unimplemented!()
}

fn parse_weight_names_line(
    line: &str,
    event: &mut Event,
) -> Result<(), ParseError> {
    let (mut rest, (_, _, nnames)) = tuple((char('N'), whitespace, int))(line)?;
    let nnames = nnames.parse()?;
    let mut weight_names = Vec::with_capacity(nnames);
    for _ in 0..nnames {
        let (rem, (_, name)) = tuple((whitespace, string))(rest)?;
        weight_names.push(name.to_owned());
        rest = rem;
    }
    event.weight_names = weight_names;
    Ok(())
}

fn parse_xs_info_line(line: &str, event: &mut Event) -> Result<(), ParseError> {
    let (_rest, parsed) = tuple((
        char('C'),
        whitespace,
        recognize_float,
        whitespace,
        recognize_float,
    ))(line)?;
    let (_, _, xs, _, xs_err) = parsed;
    let xs = CrossSection {
        cross_section: xs.parse()?,
        cross_section_error: xs_err.parse()?,
    };
    event.xs = xs;
    Ok(())
}

impl<T: BufRead> Iterator for Reader<T> {
    type Item = Result<Event, LineParseError>;

    fn next(&mut self) -> Option<Self::Item> {
        if let Err(err) = self.skip_headers() {
            return Some(Err(LineParseError {
                err: err.into(),
                line: self.line.clone(),
                line_nr: self.line_nr,
            }));
        }
        if self.line.is_empty() {
            return None;
        }
        Some(self.parse_event())
    }
}

#[derive(Debug)]
pub struct LineParseError {
    err: ParseError,
    line: String,
    line_nr: usize,
}

#[derive(Debug)]
pub enum ParseError {
    Io(io::Error),
    Parse(String),
    ConvertInt(ParseIntError),
    ConvertFloat(ParseFloatError),
    BadPrefix,
    NoVertex,
}

impl From<io::Error> for ParseError {
    fn from(err: io::Error) -> Self {
        ParseError::Io(err)
    }
}

impl From<ParseIntError> for ParseError {
    fn from(err: ParseIntError) -> Self {
        ParseError::ConvertInt(err)
    }
}

impl From<ParseFloatError> for ParseError {
    fn from(err: ParseFloatError) -> Self {
        ParseError::ConvertFloat(err)
    }
}

impl<T: Display> From<nom::Err<T>> for ParseError {
    fn from(err: nom::Err<T>) -> Self {
        match err {
            nom::Err::Failure(err) => ParseError::Parse(err.to_string()),
            nom::Err::Error(err) => ParseError::Parse(err.to_string()),
            _ => unreachable!(),
        }
    }
}

impl Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ParseError::Io(err) => write!(f, "I/O error: {}", err),
            ParseError::Parse(err) => write!(f, "Parsing error: {}", err),
            ParseError::ConvertInt(err) => {
                write!(f, "Integer conversion error: {}", err)
            }
            ParseError::ConvertFloat(err) => {
                write!(f, "Float conversion error: {}", err)
            }
            ParseError::BadPrefix => write!(f, "Unrecognized prefix"),
            ParseError::NoVertex => {
                write!(f, "Tried to add particle without vertex")
            }
        }
    }
}

impl Display for LineParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}\n in line {}:\n{}", self.err, self.line_nr, self.line)
    }
}

impl std::error::Error for LineParseError {}
