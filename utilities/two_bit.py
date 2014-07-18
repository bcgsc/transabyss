#! /usr/bin/env python

# import standard modules
import struct, sys

# import custom modules
#from error import MyError
#from messages import ErrMsg
#from files_paths import OpenFile

# CONSTANTS
TWO_BIT_SIG = 0x1A412743 # use to check whether to byte-swap
TWO_BIT_VER = 0          # only this version is currently supported
TWO_BIT_DEFAULT_TYPE = 'long'   # default field size is 4 bytes (32 bits)
ELEMENT_LENGTHS = {
  'char':  1,
  'short': 2,
  'long':  4,
}
# use format characters for unsigned values
FORMAT_CHARS = {
  'char':  'B',
  'short': 'H',
  'long':  'L',
}
TWO_BIT_BASES = {
  '00': 'T',
  '01': 'C',
  '10': 'A',
  '11': 'G',
}
BITS_PER_BASE = 2
BITS_PER_BYTE = 8

class TwoBitFileCls:
  def __init__(self, path, fail_msg="cannot open 2bit file", debug=False):
    self.debug = debug
    self.file = open(path, "rb")
    self.ReadHeader()
    self.ReadSequenceIndices()
    self.ReadSequenceHeaders()
  # end def

  def __del__(self):
    self.file.close()
  # end def

  def ResetFile(self):
    self.file.seek(0)
  # end def

  # field size in bytes
  def ReadField(self, field_type=TWO_BIT_DEFAULT_TYPE):
    if (field_type not in ELEMENT_LENGTHS):
      raise TwoBitError("invalid number type: %s. Valid types: %s" %
        (field_type, ", ".join(sorted(ELEMENT_LENGTHS.keys()))))
    # end if
    # read the number of bytes needed for one element of the specified type
    num_bytes = ELEMENT_LENGTHS[field_type]
    # bytes_str is a string of chars (each char is one byte)
    bytes_str = self.file.read(num_bytes)
    # convert to a tuple of numbers
    bytes_tup = ConvertBytesToNumbers(bytes_str, self.byte_swap, field_type)
    # return the first value of the tuple
    return bytes_tup[0]
    # if the field size was 1 byte (8 bits, 1 single char)
    #if (1 == field_type):
    #  # convert to a tuple of chars
    #  bytes_tup = ConvertBytesToNumbers(bytes_str, self.byte_swap, "char")
    #  # return the value as a single long
    #  return bytes_tup[0]
    ## if the field size was 4 bytes (32 bits, 1 single long)
    #if (TWO_BIT_FIELD_SIZE == field_type):
    #  # convert to a tuple of longs
    #  bytes_tup = ConvertBytesToNumbers(bytes_str, self.byte_swap, "long")
    #  # return the value as a single long
    #  return bytes_tup[0]
    #else:
    #  raise TwoBitError \
    #    ("invalid field size: %i byte(s). Valid field sizes: " % field_type +
    #     "1 byte (8 bits), %i bytes (%i bits)" %
    #     (TWO_BIT_FIELD_SIZE, TWO_BIT_FIELD_SIZE*8))
    # end if
  # end def

  def ReadChars(self, num_chars):
    # just read the specified number of chars (1 byte per char) and return them
    return self.file.read(num_chars)
  # end def

  def ReadArray(self, array_size, element_type=TWO_BIT_DEFAULT_TYPE):
    if (element_type not in ELEMENT_LENGTHS):
      pass
      #raise TwoBitError("invalid number type: %s. Valid types: %s" %
        #(field_type, ", ".join(sorted(ELEMENT_LENGTHS.keys()))))
    # end if
    # read enough bytes for the specified number of elements 
    # of the specified type
    num_bytes = array_size * ELEMENT_LENGTHS[element_type]
    # bytes_str is a string of chars (each char is one byte)
    bytes_str = self.file.read(num_bytes)
    # convert to a tuple of numbers
    bytes_tup = ConvertBytesToNumbers(bytes_str, self.byte_swap, element_type)
    # return the tuple
    return bytes_tup
  # end def

  def ReadHeader(self):
    # read the signature and check whether we need to byte-swap
    self.ReadSignature()
    # read the version number and check that it is valid
    self.ReadVersion()
    # read the number of sequences
    self.ReadSequenceCount()
    # skip the reserved field
    self.ReadField()
  # end def

  def ReadSignature(self):
    # assume we do not need to byte-swap
    self.byte_swap = False
    # get the signature
    sig = self.ReadField()
    if (self.debug):
      ErrMsg("signature: %i" % TWO_BIT_SIG)
      ErrMsg("sig read:  %i" % sig)
    # end if
    # check the signature (to see if we need to byte-swap)
    if (TWO_BIT_SIG != sig):
      self.ResetFile()
      self.byte_swap = True
      sig = self.ReadField()
      if (self.debug):
        ErrMsg("enabling byte-swapping")
        ErrMsg("sig read:  %i" % sig)
      # end if
      if (TWO_BIT_SIG != sig):
        pass
        #raise TwoBitError("invalid signature: %i. Should be %i." %
          #(sig, TWO_BIT_SIG))
      # end if
    # end if
  # end def

  def ReadVersion(self):
    version = self.ReadField()
    if (TWO_BIT_VER != version):
      pass
      #raise TwoBitError("invalid two-bit version: %i. Should be %i." %
        #(version, TWO_BIT_VER))
    # end if
  # end def

  def ReadSequenceCount(self):
    self.num_seqs = self.ReadField()
    if (self.debug):
      ErrMsg("Number of sequences in file: %i" % self.num_seqs)
    # end if
  # end def

  def ReadSequenceIndices(self):
    self.chrom_sequences = dict()
    #ErrMsg("ReadSequenceIndices TEMPORARILY DISABLED FOR TESTING")
    for i in xrange(self.num_seqs):
      self.ReadSequenceIndex()
    # end for
  # end def

  def ReadSequenceIndex(self):
    name_length = self.ReadField('char') # a single byte
    sequence_name = self.ReadChars(name_length)
    header_start = self.ReadField()
    #if (self.debug):
      #ErrMsg("Name length: %i" % name_length)
      #ErrMsg("Name: %s" % sequence_name)
      #ErrMsg("Offset: %i" % header_start)
    # end if
    new_sequence = \
      TwoBitSequenceCls(sequence_name, header_start, debug=self.debug)
    self.chrom_sequences[sequence_name] = new_sequence
  # end def

  def ReadSequenceHeaders(self):
    #ErrMsg("ReadSequenceHeaders TEMPORARILY DISABLED FOR TESTING")
    for sequence in self.chrom_sequences.itervalues():
      self.ReadSequenceHeader(sequence)
    # end for
  # end def

  def ReadSequenceHeader(self, sequence):
    #if (self.debug):
      #ErrMsg("Reading sequence header for %s..." % sequence.name)
    # end if
    # seek to the start of the header for the current sequence
    self.file.seek(sequence.header_start)
    # get the number of bases in the sequence
    sequence.num_bases = self.ReadField()
    # read the N-block information for the sequence
    self.ReadNBlockInfo(sequence)
    # read the mask-block information for the sequence
    self.ReadMaskBlockInfo(sequence)
    # skip the reserved field
    self.ReadField()
    # store the position in bytes from the beginning of the file of
    # the first base in the sequence
    sequence.seq_start = self.file.tell()
    # calculate the position in bytes from the beginning of the file of
    # the first base in the sequence
    #sequence.CalculateFirstSequenceBytePosition()
    # mark that the header for the sequence has been succesfully read
    sequence.header_read = True
  # end def

  def ReadNBlockInfo(self, sequence):
    # get the number of N-blocks in the sequence
    num_blocks = self.ReadField()
    #if (self.debug):
      #ErrMsg("Number of N-blocks: %i" % num_blocks)
    # end if
    # read the block starts
    sequence.n_block_starts = self.ReadArray(num_blocks)
    # read the block sizes
    sequence.n_block_sizes = self.ReadArray(num_blocks)
  # end def

  def ReadMaskBlockInfo(self, sequence):
    # get the number of mask-blocks in the sequence
    num_blocks = self.ReadField()
    #if (self.debug):
      #ErrMsg("Number of mask-blocks: %i" % num_blocks)
    # end if
    # read the block starts
    sequence.mask_block_starts = self.ReadArray(num_blocks)
    # read the block sizes
    sequence.mask_block_sizes = self.ReadArray(num_blocks)
  # end def

  def GetSequence(self, chromosome, start, end, apply_mask=False):
    # check the query region, and select the appropriate chromosome
    chrom_seq = self.CheckQueryRegion(chromosome, start, end)
    if not chrom_seq:
      return ''
    # get the starting byte position (convert 1-based to 0-based)
    (start_byte, extra_start_bases) = chrom_seq.GetStartByte(start-1)
    # get the ending byte position (convert 1-based to 0-based)
    (end_byte, extra_end_bases) = chrom_seq.GetEndByte(end-1)
    #if (self.debug):
      #ErrMsg("Start byte: %i" % start_byte)
      #ErrMsg("Extra start bases: %i" % extra_start_bases)
      #ErrMsg("End byte: %i" % end_byte)
      #ErrMsg("Extra end bases: %i" % extra_end_bases)
    # end if
    # get the sequence fragment
    seq_frag = self.GetSequenceFragment(start_byte, end_byte)
    # trim the sequence fragment as required
    seq_frag = self.TrimSequence(seq_frag, extra_start_bases, extra_end_bases)
    seq_frag = chrom_seq.ApplyNBlocks(seq_frag, start-1, end-1)
    if (apply_mask):
      seq_frag = chrom_seq.ApplyMaskBlocks(seq_frag, start-1, end-1)
    # end if
    return seq_frag
  # end def

  def CheckQueryRegion(self, chromosome, start, end):
    errors = list()
    # check that chromosome is a valid sequence identifier
    chrom_valid = True
    if (chromosome not in self.chrom_sequences):
      # if it does not already start with "chr", try adding "chr"
      if (not chromosome.startswith("chr")):
        if ("chr%s" % chromosome in self.chrom_sequences):
          chromosome = "chr%s" % chromosome
        else:
          chrom_valid = False
        # end if
      else:
        chrom_valid = False
      # end if
      if (not chrom_valid):
        errors.append("invalid chromosome: %s." % chromosome)
      # end if
    # end if
    # check that start and end are valid base coordinates
    if (1 > start):
      errors.append("invalid start coordinate: "
        "%i. Must be 1 or greater." % start)
    # end if
    if (end < start):
      errors.append("invalid coordinates: "
        "start (%i) must be less than or equal to end (%i)." % (start, end))
    # end if
    chrom_seq = None
    if (chrom_valid):
      chrom_seq = self.chrom_sequences[chromosome]
      # check that sequence has had its header read
      if (not chrom_seq.header_read):
        self.ReadSequenceHeader(chrom_seq)
      # end if
      if (chrom_seq.num_bases < end):
        errors.append("invalid end coordinate: %i. Must be %i or less." %
          (end, chrom_seq.num_bases))
      # end if
    # end if
    if (0 < len(errors)):
      pass
      #sys.stderr.write("Invalid sequence region query:\n"
        #"%s" % "\n".join(errors))
      #raise TwoBitError("Invalid sequence region query:\n"
        #"%s" % "\n".join(errors))
    # end if
    return chrom_seq
  # end def

  def GetSequenceFragment(self, start_byte, end_byte):
    # move to the first byte
    self.file.seek(start_byte)
    # read the correct number of bytes
    num_bytes = end_byte - start_byte + 1
    bytes = self.ReadChars(num_bytes)
    # convert the bytes to a bit string
    bits = BytesToBitString(bytes)
    # convert the bit string to bases
    bases = BitStringToBaseString(bits)
    #if (self.debug):
      #ErrMsg("Full sequence fragment:\n%s" % bases)
    # end if
    # return the bases
    return bases
  # end def

  def TrimSequence(self, seq, start_trim, end_trim):
    if (0 == end_trim):
      trimmed_seq = seq[start_trim:]
    else:
      trimmed_seq = seq[start_trim:-end_trim]
    # end if
    #if (self.debug):
      #ErrMsg("Trimmed sequence:\n%s" % trimmed_seq)
    # end if
    return trimmed_seq
  # end def
# end class

class TwoBitSequenceCls:
  def __init__(self, name, header_start, debug=False):
    self.name   = name
    self.header_start = header_start
    self.debug  = debug
    self.header_read = False
    # num_bases
    # n_block_starts
    # n_block_sizes
    # mask_block_starts
    # mask_block_sizes
    # seq_start
  # end def

  def CalculateFirstSequenceBytePosition(self):
    # start with the initial header_start
    self.seq_start = self.header_start
    # add an offset for the number of bases field
    self.seq_start += ELEMENT_LENGTHS[TWO_BIT_DEFAULT_TYPE]
    # add an offset for the N-block count field
    self.seq_start += ELEMENT_LENGTHS[TWO_BIT_DEFAULT_TYPE]
    # add an offset for the N-block starts array
    self.seq_start += (len(self.n_block_starts) *
                       ELEMENT_LENGTHS[TWO_BIT_DEFAULT_TYPE])
    # add an offset for the N-block sizes array
    self.seq_start += (len(self.n_block_sizes) *
                       ELEMENT_LENGTHS[TWO_BIT_DEFAULT_TYPE])
    # add an offset for the mask-block count field
    self.seq_start += ELEMENT_LENGTHS[TWO_BIT_DEFAULT_TYPE]
    # add an offset for the mask-block starts array
    self.seq_start += (len(self.mask_block_starts) *
                       ELEMENT_LENGTHS[TWO_BIT_DEFAULT_TYPE])
    # add an offset for the mask-block sizes array
    self.seq_start += (len(self.mask_block_sizes) *
                       ELEMENT_LENGTHS[TWO_BIT_DEFAULT_TYPE])
    # add an offset for the reserved field
    self.seq_start += ELEMENT_LENGTHS[TWO_BIT_DEFAULT_TYPE]
  # end def

  def GetStartByte(self, start_base):
    start_bit  = start_base * BITS_PER_BASE
    start_byte = start_bit / BITS_PER_BYTE
    extra_start_bits  = start_bit - (start_byte * BITS_PER_BYTE)
    extra_start_bases = extra_start_bits / BITS_PER_BASE
    if (BITS_PER_BYTE / BITS_PER_BASE <= extra_start_bases):
      pass
      #raise TwoBitError("error getting start byte for base in "
        #"position %i, too many extra bases: %i" %
        #(start_base+1, extra_start_bases))
    # end if
    return (self.seq_start+start_byte, extra_start_bases)
  # end def

  def GetEndByte(self, end_base):
    end_bit  = (end_base*BITS_PER_BASE) + (BITS_PER_BASE-1)
    end_byte = end_bit / BITS_PER_BYTE
    extra_end_bits  = BITS_PER_BYTE*(end_byte+1) - (end_bit+1)
    extra_end_bases = extra_end_bits / BITS_PER_BASE
    if (BITS_PER_BYTE / BITS_PER_BASE <= extra_end_bases):
      pass
      #raise TwoBitError("error getting end byte for base in "
        #"position %i, too many extra bases: %i" %
        #(end_base+1, extra_end_bases))
    # end if
    return (self.seq_start+end_byte, extra_end_bases)
  # end def

  def ApplyNBlocks(self, seq, start, end):
    #if (self.debug):
      #ErrMsg("Applying N blocks...")
    # end if
    return self.ApplyBlocks(seq, start, end,
      self.n_block_starts, self.n_block_sizes, self.ApplyNBlock)
  # end def

  def ApplyMaskBlocks(self, seq, start, end):
    #if (self.debug):
      #ErrMsg("Applying mask blocks...")
    # end if
    return self.ApplyBlocks(seq, start, end,
      self.mask_block_starts, self.mask_block_sizes, self.ApplyMaskBlock)
  # end def

  def ApplyBlocks(self, seq, start, end, block_starts, block_sizes,
      ApplyBlock):
    #if (self.debug):
    #  ErrMsg("Seq start: %i, end: %i" % (start, end))
    # end if
    new_seq = seq
    for (block_start, block_size) in zip(block_starts, block_sizes):
      block_end = block_start + block_size - 1
      #if (self.debug):
      #  ErrMsg("Block start: %i, end: %i (size: %i)" %
      #    (block_start, block_end, block_size))
      # end if
      # skip blocks that end before or start after the sequence
      if (block_end < start or end < block_start):
        #if (self.debug):
        #  ErrMsg("--Skipping block")
        # end if
        continue
      # end if
      # adjust start position if block starts before sequence
      if (block_start < start):
        block_start = start
      # end if
      # adjust end position if block ends after sequence
      if (end < block_end):
        block_end = end
      # end if
      # convert the genomic coordinates to positions in the sequence string
      string_start = block_start - start
      string_end   = block_end   - start
      # apply the block
      new_seq = ApplyBlock(new_seq, string_start, string_end)
    # end for
    return new_seq
  # end def

  def ApplyNBlock(self, seq, start, end):
    if (self.debug):
      ErrMsg("Applying N block to %i-%i" % (start, end))
    # end if
    block_size = end - start + 1
    new_seq = seq[:start]
    new_seq += "N"*block_size
    new_seq += seq[end+1:]
    if (len(seq) != len(new_seq)):
      pass
      #raise TwoBitError("error applying N-block.\nOld seq:%s\nNew seq:%s" %
        #(seq, new_seq))
    # end if
    return new_seq
  # end def

  def ApplyMaskBlock(self, seq, start, end):
    if (self.debug):
      ErrMsg("Applying mask block to %i-%i" % (start, end))
    # end if
    block_size = end - start + 1
    new_seq = seq[:start]
    new_seq += seq[start:end+1].lower()
    new_seq += seq[end+1:]
    if (len(seq) != len(new_seq)):
      pass
      #raise TwoBitError("error applying mask-block.\nOld seq:%s\nNew seq:%s" %
        #(seq, new_seq))
    # end if
    return new_seq
  # end def
# end class

def ConvertBytesToNumbers(bytes_str, byte_swap, element_type):
  if (byte_swap):
    endian = "<"
  else:
    endian = ">"
  # end if
  if (element_type not in ELEMENT_LENGTHS):
    pass
    #raise TwoBitError("invalid number type: %s. Valid types: %s" %
      #(element_type, ", ".join(sorted(ELEMENT_LENGTHS.keys()))))
  # end if
  value_len = ELEMENT_LENGTHS[element_type]
  num_values = len(bytes_str) / value_len
  if (num_values * value_len != len(bytes_str)):
    pass
    #raise TwoBitError("invalid byte string, the length of the byte string "
      #"(%i) " % len(bytes_str) + "must be a multiple of the byte_length of "
      #"the number type being converted to (%i)" % value_len)
  # end if
  format_char = FORMAT_CHARS[element_type]
  pattern = "%s%i%s" % (endian, num_values, format_char)
  int_values_tuple = struct.unpack(pattern, bytes_str)
  return int_values_tuple
# end def

# need to test whether this needs to be changed if the file
# needs to be byte-swapped (does not look like it though)
def BytesToBitString(bytes_str):
  bit_string = ""
  # bytes_str will be a string (of chars)
  for byte_char in bytes_str:
    # convert to an integer
    byte_int = ord(byte_char)
    # get each of the 8 bits of the byte
    for i in reversed(xrange(8)):
      bit_string += "%i" % ((byte_int >> i) & 1)
    # end for
  # end for
  return bit_string
# end def

def BitStringToBaseString(bit_string):
  base_string = ""
  for i in xrange(0, len(bit_string), 2):
    two_bits = bit_string[i:i+2]
    base_string += TWO_BIT_BASES[two_bits]
  # end for
  return base_string
# end def

def TestGetSeq(file=None,bytes_str=None):
  if (None != file):
    bytes_str = file.read(13)
  elif(None == bytes_str):
    ErrMsg("use one of file or bytes_str")
    return
  # end if
  print bytes_str
  print "No reversing..."
  bit_string = ""
  # bytes_str will be a string (of chars)
  for byte_char in bytes_str:
    # convert to an integer
    byte_int = ord(byte_char)
    # get each of the 8 bits of the byte
    for i in xrange(8):
      bit_string += "%i" % ((byte_int >> i) & 1)
    # end for
  # end for
  print BitStringToBaseString(bit_string)
  print "Reversing bits..."
  bit_string = ""
  # bytes_str will be a string (of chars)
  for byte_char in bytes_str:
    # convert to an integer
    byte_int = ord(byte_char)
    # get each of the 8 bits of the byte
    for i in reversed(xrange(8)):
      bit_string += "%i" % ((byte_int >> i) & 1)
    # end for
  # end for
  print BitStringToBaseString(bit_string)
  print "Reversing bytes..."
  bit_string = ""
  # bytes_str will be a string (of chars)
  for byte_char in reversed(bytes_str):
    # convert to an integer
    byte_int = ord(byte_char)
    # get each of the 8 bits of the byte
    for i in xrange(8):
      bit_string += "%i" % ((byte_int >> i) & 1)
    # end for
  # end for
  print BitStringToBaseString(bit_string)
  print "Reversing both..."
  bit_string = ""
  # bytes_str will be a string (of chars)
  for byte_char in reversed(bytes_str):
    # convert to an integer
    byte_int = ord(byte_char)
    # get each of the 8 bits of the byte
    for i in reversed(xrange(8)):
      bit_string += "%i" % ((byte_int >> i) & 1)
    # end for
  # end for
  print BitStringToBaseString(bit_string)
# end def

#### EXCEPTION CLASSES ####
class TwoBitError():
  """Exception raised for errors encountered by two bit classes"""
  pass
# end class
