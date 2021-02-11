use std::collections::HashMap;
use std::default::Default;
use std::fmt::{self, Display};
use std::io::{self, BufRead};
use std::iter::Iterator;
use std::num::{ParseIntError, ParseFloatError};

use crate::event::*;

use nom::{
    character::{
        complete::{char, space1}
    },
    bytes::complete::{take_until, take_while1},
    combinator::{opt, recognize},
    number::complete::{double, recognize_float},
    sequence::{delimited, preceded, tuple},
    IResult
};

const BUF_SIZE: usize = 256;

/// Reader for the HepMC2 format
#[derive(Debug, Eq, PartialEq, Ord, PartialOrd, Hash, Default)]
pub struct Reader<T> {
    stream: T,
    line: String,
    line_nr: usize,
}

impl<T: BufRead + Default> Reader<T> {
    pub fn new() -> Self {
        Self::default()
    }
}

impl<T: BufRead> From<T> for Reader<T> {
    fn from(stream: T) -> Self {
        Self{stream, line: String::with_capacity(BUF_SIZE), line_nr: 0}
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
                break
            };
            self.line_nr += 1;
            match self.line.as_bytes().first() {
                Some(b'E') => break,
                Some(b'V') => parse_vertex_line(&self.line, &mut event)?,
                Some(b'P') => parse_particle_line(&self.line, &mut event)?,
                Some(b'U') => parse_units_line(&self.line, &mut event)?,
                Some(b'F') => parse_pdf_info_line(&self.line, &mut event)?,
                Some(b'H') => parse_heavy_ion_line(&self.line, &mut event)?,
                Some(b'N') => parse_weight_names_line(&self.line, &mut event)?,
                Some(b'C') => parse_xs_info_line(&self.line, &mut event)?,
                _ => return Err(ParseError::BadPrefix)
            };
        }
        Ok(event)
    }

    fn parse_event(&mut self) -> Result<Event, LineParseError> {
        self.parse_event_inner().map_err(
            |err| LineParseError{
                err,
                line: self.line.clone(),
                line_nr: self.line_nr,
            }
        )
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
        take_while1(|c: char| c.is_ascii_digit())
    )));
    int_parser(line)
}

fn string(line: &str) -> IResult<&str, &str> {
    delimited(char('"'), take_until("\""), char('"'))(line)
}

fn parse_event_line(line: &str) -> Result<Event, ParseError> {
    let (rest, parsed) = tuple((
        char('E'), whitespace,
        int, whitespace, // event number
        int, whitespace, // mpi
        recognize_float, whitespace, //event scale
        recognize_float, whitespace, // alpha_qcd
        recognize_float, whitespace, // alpha_qed
        int, whitespace, // signal_process_id
        int, whitespace, // signal_process_vertex
        int, whitespace, // num_vertices
        int, whitespace, // beam1
        // have to stop here, tuple can't handle more entries
    ))(line)?;
    let (
        _, _,
        event_number, _,
        mpi, _,
        event_scale, _,
        alpha_qcd, _,
        alpha_qed, _,
        signal_process_id, _,
        signal_process_vertex, _,
        num_vertices, _,
        _beam1, _,
    ) = parsed;
    let (mut rest, parsed) = tuple((
        int, whitespace, // beam2
        int, // random states
    ))(rest)?;
    let (_beam2, _, nrandom_states) = parsed;
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
        let (rem, weight) = preceded(whitespace, double)(rest)?;
        rest = rem;
        weights.push(weight);
    }
    let event = Event{
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
        weight_names:Default::default(),
        xs:          Default::default(),
        energy_unit: Default::default(),
        length_unit: Default::default(),
        pdf_info:    Default::default(),
    };
    Ok(event)
}

fn parse_vertex_line(line: &str, event: &mut Event) -> Result<(), ParseError> {
    let (mut rest, parsed) = tuple((
        char('V'), whitespace,
        int, whitespace, // barcode
        int, whitespace, // status
        recognize_float, whitespace, // x
        recognize_float, whitespace, // y
        recognize_float, whitespace, // z
        recognize_float, whitespace, // t
        int, whitespace, // num_orphans_in
        int, whitespace, // num_particles_out
        int, // num_weights
    ))(line)?;
    let (
        _, _,
        barcode, _,
        status, _,
        x, _,
        y, _,
        z, _,
        t, _,
        _num_orphans_in, _,
        num_particles_out, _,
        num_weights,
    ) = parsed;
    let num_weights = num_weights.parse()?;
    let mut weights = Vec::with_capacity(num_weights);
    for _ in 0..num_weights {
        let (rem, weight) = preceded(whitespace, double)(rest)?;
        rest = rem;
        weights.push(weight);
    }
    let vertex = Vertex{
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

fn parse_particle_line(line: &str, event: &mut Event) -> Result<(), ParseError> {
    let (rest, parsed) = tuple((
        char('P'), whitespace,
        int, whitespace, // barcode
        int, whitespace, // id
        recognize_float, whitespace, // px
        recognize_float, whitespace, // py
        recognize_float, whitespace, // pz
        recognize_float, whitespace, // E
        recognize_float, whitespace, // m
    ))(line)?;
    let (
        _, _,
        _barcode, _,
        id, _,
        px, _,
        py, _,
        pz, _,
        e, _,
        m, _,
    ) = parsed;
    let (mut rest, parsed) = tuple((
        int, whitespace, // status
        recognize_float, whitespace, // theta
        recognize_float, whitespace, // phi
        int, whitespace, // end_vtx_code
        int, // flowsize
    ))(rest)?;
    let (
        status, _,
        theta, _,
        phi, _,
        end_vtx_code, _,
        flowsize
    ) = parsed;
    let mut flows = HashMap::new();
    for _ in 0..flowsize.parse()? {
        let (rem, flowidx) = preceded(whitespace, int)(rest)?;
        let flowidx = flowidx.parse()?;
        let (rem, flowval) = preceded(whitespace, int)(rem)?;
        let flowval = flowval.parse()?;
        rest = rem;
        flows.insert(flowidx, flowval);
    }
    let particle = Particle{
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
        return Err(ParseError::NoVertex)
    }
    Ok(())
}

fn parse_units_line(line: &str, event: &mut Event) -> Result<(), ParseError> {
    let (_rest, parsed) = tuple((
        char('U'), whitespace,
        non_whitespace, whitespace,
        non_whitespace
    ))(line)?;
    let (_, _, energy, _, length) = parsed;
    event.energy_unit = energy.to_owned();
    event.length_unit = length.to_owned();
    Ok(())
}

fn parse_pdf_info_line(line: &str, event: &mut Event) -> Result<(), ParseError> {
    let (_rest, parsed) = tuple((
        char('F'), whitespace,
        int, whitespace, // id0
        int, whitespace, // id1
        recognize_float, whitespace, // x0
        recognize_float, whitespace, // x1
        recognize_float, whitespace, // scale
        recognize_float, whitespace, // xf0
        recognize_float, whitespace, // xf1
        opt(int), whitespace, // pdf_id0
        opt(int), whitespace, // pdf_id1
    ))(line)?;
    let (
        _, _,
        id0, _,
        id1, _,
        x0, _,
        x1, _,
        scale, _,
        xf0, _,
        xf1, _,
        pdf_id0, _,
        pdf_id1, _,
    ) = parsed;
    let pdf_info = PdfInfo{
        parton_id: [id0.parse()?, id1.parse()?],
        x: [x0.parse()?, x1.parse()?],
        scale: scale.parse()?,
        xf: [xf0.parse()?, xf1.parse()?],
        pdf_id: [
            pdf_id0.map(|id| id.parse()).transpose()?.unwrap_or(0),
            pdf_id1.map(|id| id.parse()).transpose()?.unwrap_or(0)
        ]
    };
    event.pdf_info = pdf_info;
    Ok(())
}

fn parse_heavy_ion_line(_line: &str, _event: &mut Event) -> Result<(), ParseError> {
    unimplemented!()
}

fn parse_weight_names_line(line: &str, event: &mut Event) -> Result<(), ParseError> {
    let (mut rest, (_, _, nnames)) = tuple((
        char('N'), whitespace,
        int
    ))(line)?;
    let nnames = nnames.parse()?;
    let mut weight_names = Vec::with_capacity(nnames);
    for _ in 0..nnames {
        let (rem, (_, name)) = tuple((
            whitespace, string
        ))(rest)?;
        weight_names.push(name.to_owned());
        rest = rem;
    }
    event.weight_names = weight_names;
    Ok(())
}

fn parse_xs_info_line(line: &str, event: &mut Event) -> Result<(), ParseError> {
    let (_rest, parsed) = tuple((
        char('C'), whitespace,
        recognize_float, whitespace,
        recognize_float,
    ))(line)?;
    let (_, _, xs, _, xs_err) = parsed;
    let xs = CrossSection{
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
            return Some(Err(LineParseError{
                err: err.into(),
                line: self.line.clone(),
                line_nr: self.line_nr
            }))
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
    fn from(err: io::Error)  -> Self {
        ParseError::Io(err)
    }
}

impl From<ParseIntError> for ParseError {
    fn from(err: ParseIntError)  -> Self {
        ParseError::ConvertInt(err)
    }
}

impl From<ParseFloatError> for ParseError {
    fn from(err: ParseFloatError)  -> Self {
        ParseError::ConvertFloat(err)
    }
}

impl<T: Display> From<nom::Err<T>> for ParseError {
    fn from(err: nom::Err<T>)  -> Self {
        match err {
            nom::Err::Failure(err) => ParseError::Parse(err.to_string()),
            nom::Err::Error(err) => ParseError::Parse(err.to_string()),
            _ => unreachable!()
        }
    }
}

impl Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ParseError::Io(err) => write!(f, "I/O error: {}", err),
            ParseError::Parse(err) => write!(f, "Parsing error: {}", err),
            ParseError::ConvertInt(err) => write!(f, "Integer conversion error: {}", err),
            ParseError::ConvertFloat(err) => write!(f, "Float conversion error: {}", err),
            ParseError::BadPrefix => write!(f, "Unrecognized prefix"),
            ParseError::NoVertex => write!(f, "Tried to add particle without vertex"),
        }
    }
}

impl Display for LineParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}\n in line {}:\n{}", self.err, self.line_nr, self.line)
    }
}

impl std::error::Error for LineParseError { }

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tst_parse() {
        let line = r#"
HepMC::Version 2.06.09
HepMC::IO_GenEvent-START_EVENT_LISTING
E 0 -1 -1.0000000000000000e+00 -1.0000000000000000e+00 -1.0000000000000000e+00 0 0 23 1 2 0 16 5.5606031127834702e-09 5.5606031127834702e-09 0 -1.0000000000000000e+00 2.6206847299999998e+00 5.8670509300000004e+00 6.1711295699999997e+00 -1.0000000000000000e+00 5.9781127600000001e+00 3.9987750200000001e+00 2.8889888799999999e+00 5.5606031100000000e-09 -1.0000000000000000e+00 8.6938853699999996e-01 3.5507575000000002e-01 5.6137609399999999e-01
N 16 "0" "Weight" "eventNumber" "phi" "phi1" "phi2" "phi3" "t" "t1" "t2" "t3" "unlops_weight" "z" "z1" "z2" "z3"
U GEV MM
C 5.5606031127834701e+00 5.3451183000000005e+04
V -1 0 0 0 0 0 0 1 0
P 3 -1 0 0 5.2051533588697652e+01 5.2051533588697652e+01 0 21 0 0 -3 1 2 501
V -2 0 0 0 0 0 0 2 0
P 4 2 0 0 -2.9818545945620773e+01 2.9818545945620773e+01 0 21 0 0 -3 1 1 501
P 9 21 5.7076480059072443e+00 3.2748659206013784e+00 -2.4249351062062976e+00 7.0130095413142710e+00 0 43 0 0 -8 2 1 502 2 501
V -3 0 0 0 0 0 0 1 0
P 5 24 -4.6339895476421589e-08 -1.0408832373798305e-07 2.2232987710669342e+01 8.1870080240182716e+01 7.8793428000000006e+01 22 0 0 -6 0
V -4 0 0 0 0 0 1 1 0
P 1 2212 0 0 6.9999999371178146e+03 7.0000000000000000e+03 9.3827000000000005e-01 4 0 0 -4 0
P 6 -1 0 0 5.2051533588697652e+01 5.2051533588697652e+01 0 42 0 0 -1 1 2 501
V -5 0 0 0 0 0 0 2 0
P 7 2 0 0 -3.6129848231502301e+01 3.6129848231502301e+01 0 41 0 0 -2 1 1 502
P 12 -2 2.5237885213344935e+00 3.3116569526639856e+00 -6.9143141497758080e+01 6.9268395365069210e+01 0 1 0 0 0 1 2 503
V -6 0 0 0 0 0 0 1 0
P 8 24 -5.7076481044982819e+00 -3.2748660546697370e+00 1.8346620499347921e+01 8.1168372992804763e+01 7.8793428000000006e+01 44 0 0 -23 0
V -7 0 0 0 0 0 0 2 0
P 10 21 0 0 -1.0175250837889556e+02 1.0175250837889556e+02 0 41 0 0 -5 2 1 502 2 503
P 15 21 2.9123097878176076e+00 2.8941914776608224e-01 -8.0591874171438874e+01 8.0644996709085106e+01 0 43 0 0 -18 2 1 503 2 504
V -8 0 0 0 0 0 0 1 0
P 11 21 3.1838594845727508e+00 -3.6791032062607254e-02 1.0955462441585269e+00 3.3672743236383127e+00 0 48 0 0 -10 2 1 502 2 501
V -9 0 0 0 0 0 0 2 0
P 13 21 0 0 -1.8125550661114170e+02 1.8125550661114170e+02 0 41 0 0 -7 2 1 502 2 504
P 18 2 -6.0739028889597069e-01 2.5046574652089002e+00 -2.2999284645202547e+03 2.2999299085281878e+03 0 43 0 0 -17 1 1 504
V -10 0 0 0 0 0 0 1 0
P 14 21 2.7154969675514318e-01 -3.2621017982868949e-01 2.1844221833512592e+00 2.2252758467993488e+00 0 48 0 0 -12 2 1 502 2 501
V -11 0 0 0 0 0 0 2 0
P 16 2 0 0 -2.4821608453435269e+03 2.4821608453435274e+03 0 41 0 0 -9 1 1 502
P 21 21 -1.5176954196303214e+00 -1.1581464938267174e+00 -2.1424758388954547e+02 2.1425608954052109e+02 0 43 0 0 -16 2 1 505 2 502
V -12 0 0 0 0 0 0 1 0
P 17 21 8.7893998565111386e-01 -2.8308676450375896e+00 1.2075479712207198e+00 3.2007060509972689e+00 0 48 0 0 -14 2 1 502 2 501
V -13 0 0 0 0 0 1 1 0
P 2 2212 0 0 -6.9999999371178146e+03 7.0000000000000000e+03 9.3827000000000005e-01 4 0 0 -13 0
P 19 2 0 0 -2.6963868241367072e+03 2.6963868241367081e+03 0 41 0 0 -11 1 1 505
V -14 0 0 0 0 0 0 1 0
P 20 21 2.3966354052814354e+00 -1.6727211512108722e+00 1.2291530675859121e+00 3.1705953036569099e+00 0 48 0 0 -15 2 1 502 2 501
V -15 0 0 0 0 0 0 2 0
P 22 21 1.6780099262880444e+00 -1.2585715889494420e-02 1.5618589221168571e+00 2.2924395315359480e+00 0 1 0 0 0 2 1 506 2 501
P 23 21 7.0154368522484822e-01 -1.6731704747007090e+00 -2.7440810198086751e+00 3.2896266692386575e+00 0 1 0 0 0 2 1 502 2 506
V -16 0 0 0 0 0 0 1 0
P 24 21 -1.5006136258617786e+00 -1.1451114544473862e+00 -2.1183620872426775e+02 2.1184461864340338e+02 0 1 0 0 0 2 1 505 2 502
V -17 0 0 0 0 0 0 2 0
P 25 2 -1.6033704065363041e+00 4.3863302327067144e+00 -2.1547456957401910e+03 2.1547507568163073e+03 0 51 0 0 -20 1 1 507
P 26 21 2.6011093609332741e+00 -1.7221584398720926e+00 -1.8960124875750688e+02 1.8962691035295006e+02 0 51 0 0 -19 2 1 504 2 507
V -18 0 0 0 0 0 0 1 0
P 27 21 1.3071805445246669e+00 1.2990482014036053e-01 -3.6173394193995875e+01 3.6197238068015622e+01 0 1 0 0 0 2 1 503 2 504
V -19 0 0 0 0 0 0 2 0
P 28 21 3.9678447216374000e+00 -1.0579393283155296e+00 -3.5312021479758994e+02 3.5314409116683936e+02 0 51 0 0 -22 2 1 504 2 508
P 29 21 -1.9866544379963598e+00 1.0316895689562262e+00 -6.6958108247799419e+02 6.6958482449322673e+02 0 51 0 0 -21 2 1 508 2 507
V -20 0 0 0 0 0 0 1 0
P 30 2 -9.8345132924407030e-01 2.6904215521939250e+00 -1.3216456472221134e+03 1.3216487515091915e+03 0 1 0 0 0 1 1 507
V -21 0 0 0 0 0 0 2 0
P 31 21 -9.8001674572977349e-01 1.3322582893034975e+00 -6.1464412971562535e+02 6.1464635485688007e+02 0 1 0 0 0 2 1 509 2 507
P 32 21 -6.7009385001849431e-01 -3.9030080197012751e-01 -8.4887830929292363e+01 8.4891372944274849e+01 0 1 0 0 0 2 1 508 2 509
V -22 0 0 0 0 0 0 1 0
P 33 21 3.6313008793893080e+00 -9.6820724669267344e-01 -3.2316933663066641e+02 3.2319118785891123e+02 0 1 0 0 0 2 1 504 2 508
V -23 0 0 0 0 0 0 2 0
P 34 -13 6.7280635852793118e+00 2.6135691981490424e+01 3.5757112587329956e+01 4.4798699749357091e+01 1.0565837000000000e-01 1 0 0 0 0
P 35 14 -1.2435711898718612e+01 -2.9410557911101773e+01 -1.7410492121416357e+01 3.6369672571802241e+01 0 1 0 0 0 0
"#;
        let mut reader = Reader::from(line.as_bytes());
        assert!(reader.next().unwrap().is_ok());
        assert!(reader.next().is_none());
    }
}
