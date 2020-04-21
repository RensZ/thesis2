PDS_VERSION_ID               = PDS3
RECORD_TYPE                  = FIXED_LENGTH
RECORD_BYTES                 = 122
FILE_RECORDS                 = 1327
^SHADR_HEADER_TABLE          = ("GGMES_100V07_SHA.TAB",1)
^SHADR_COEFFICIENTS_TABLE    = ("GGMES_100V07_SHA.TAB",3)
INSTRUMENT_HOST_NAME         = "MESSENGER"
TARGET_NAME                  = "MERCURY"
INSTRUMENT_NAME              = "RADIO SCIENCE SUBSYSTEM"
DATA_SET_ID                  = "MESS-H-RSS/MLA-5-SDP-V1.0"
OBSERVATION_TYPE             = "GRAVITY FIELD"
ORIGINAL_PRODUCT_ID          = "HGM007A.SHA"
PRODUCT_ID                   = "HGM007A.SHA"
PRODUCT_RELEASE_DATE         = 2016-05-06
DESCRIPTION                  = "This file contains coefficient
and related data for a spherical harmonic model of the Hermean
gravity field.  Input data are from radio tracking of the MErcury
Surface, Space ENvironment, GEochemistry, and Ranging (MESSENGER)
spacecraft.  This product is a set of binary tables:
a header table, a names table, a coefficients table, and a
covariance table.  Definitions of the tables follow.  This
MESSENGER Mercury gravity model is in the form of a Spherical
Harmonics Ascii Data Record (SHADR).  It has been produced by
the MESSENGER Radio Science at NASA GSFC, and extends to degree
and order 100. A Kaula rule of 3.0E-05/L^2 was applied.
All the MESSENGER tracking data, up to the end of mission
(April 30, 2015) are included in this solution.

The Mercury orientation was chosen to follow the parameters
adopted by the MESSENGER project for the final PDS delivery, as
described in the PCK kernel 'pck00010_msgr_v21.tpc'.
The pole orientation, prime meridian, and spin period were
not adjusted during the gravity solution.

Reference radius = 2439.4 km

Pole orientation (right ascension and declination at J2000):
 RA = 281.01030 degree
DEC =  61.41556 degree

Prime meridian = 329.5988 degrees (at J2000)
Spin rate = 6.1385108 degrees/day.
  
The citation for this model is MAZARICOETAL2016."

START_ORBIT_NUMBER           = 1
STOP_ORBIT_NUMBER            = 4104
START_TIME                   = 2008-01-07
STOP_TIME                    = 2015-04-30
PRODUCT_CREATION_TIME        = 2016-01-22
PRODUCER_FULL_NAME           = "ERWAN MAZARICO"
PRODUCER_INSTITUTION_NAME    = "GODDARD SPACE FLIGHT CENTER"
PRODUCER_ID                  = "MESSENGER RS TEAM"
SOFTWARE_NAME                = "NSOLVE"

OBJECT                     = SHADR_HEADER_TABLE
ROWS                         = 1
COLUMNS                      = 8
ROW_BYTES                    = 137
ROW_SUFFIX_BYTES             = 107
INTERCHANGE_FORMAT           = ASCII
DESCRIPTION                  = "The SHADR header includes descriptive
 information about the spherical harmonic coefficients which follow in
 SHADR_COEFFICIENTS_TABLE.  The header consists of a single record of eight
 (delimited) data columns requiring 137 bytes, a pad of 105 unspecified ASCII
 characters, an ASCII carriage-return, and an ASCII line-feed."

  OBJECT                   = COLUMN
    NAME                         = "REFERENCEADIUS"
    DATA_TYPE                    = ASCII_REAL
    START_BYTE                   = 1
    BYTES                        = 23
    FORMAT                       = "E23.16"
    UNIT                         = "KILOMETER"
    DESCRIPTION                  = "The assumed reference radius of the
     spherical planet."
  END_OBJECT               = COLUMN

  OBJECT                   = COLUMN
    NAME                         = "CONSTANT"
    DATA_TYPE                    = ASCII_REAL
    START_BYTE                   = 25
    BYTES                        = 23
    FORMAT                       = "E23.16"
    UNIT                         = "N/A"
    DESCRIPTION                  = "For a gravity field model the assumed
     gravitational constant GM in km cubed per seconds squared for the
     planet.  For a topography model, set to 1."
  END_OBJECT               = COLUMN

  OBJECT                   = COLUMN
    NAME                         = "UNCERTAINTY IN CONSTANT"
    DATA_TYPE                    = ASCII_REAL
    START_BYTE                   = 49
    BYTES                        = 23
    FORMAT                       = "E23.16"
    UNIT                         = "N/A"
    DESCRIPTION                  = "For a gravity field model the
     uncertainty in the gravitational constant GM in km cubed per second
     squared for the planet (or, set to 0 if not known). For a topography
     model, set to 0."
  END_OBJECT               = COLUMN

  OBJECT                   = COLUMN
    NAME                         = "DEGREE OF FIELD"
    DATA_TYPE                    = ASCII_INTEGER
    START_BYTE                   = 73
    BYTES                        = 5
    FORMAT                       = "I5"
    UNIT                         = "N/A"
    DESCRIPTION                  = "Degree of the model field."
  END_OBJECT               = COLUMN

  OBJECT                   = COLUMN
    NAME                         = "ORDER OF FIELD"
    DATA_TYPE                    = ASCII_INTEGER
    START_BYTE                   = 79
    BYTES                        = 5
    FORMAT                       = "I5"
    UNIT                         = "N/A"
    DESCRIPTION                  = "Order of the model field."
  END_OBJECT               = COLUMN

  OBJECT                   = COLUMN
    NAME                         = "NORMALIZATION STATE"
    DATA_TYPE                    = ASCII_INTEGER
    START_BYTE                   = 85
    BYTES                        = 5
    FORMAT                       = "I5"
    UNIT                         = "N/A"
    DESCRIPTION                  = "The normalization indicator.
     For gravity field:
        0   coefficients are unnormalized
        1   coefficients are normalized
        2   other."
  END_OBJECT               = COLUMN

  OBJECT                   = COLUMN
    NAME                         = "REFERENCE LONGITUDE"
    POSITIVE_LONGITUDE_DIRECTION = "EAST"
    DATA_TYPE                    = ASCII_REAL
    START_BYTE                   = 91
    BYTES                        = 23
    FORMAT                       = "E23.16"
    UNIT                         = "DEGREE"
    DESCRIPTION                  = "The reference longitude for the
     spherical harmonic expansion; normally 0."
  END_OBJECT               = COLUMN

  OBJECT                   = COLUMN
    NAME                         = "REFERENCE LATITUDE"
    DATA_TYPE                    = ASCII_REAL
    START_BYTE                   = 115
    BYTES                        = 23
    FORMAT                       = "E23.16"
    UNIT                         = "DEGREE"
    DESCRIPTION                  = "The reference latitude for the
     spherical harmonic expansion; normally 0."
  END_OBJECT               = COLUMN

END_OBJECT                 = SHADR_HEADER_TABLE

OBJECT                     = SHADR_COEFFICIENTS_TABLE
  ROWS                       = 5150
  COLUMNS                    = 6
  ROW_BYTES                  = 107
  ROW_SUFFIX_BYTES           = 15
  INTERCHANGE_FORMAT         = ASCII
  DESCRIPTION                = "The SHADR coefficients table contains the
   coefficients for the spherical harmonic model. Each row in the table
   contains the degree index m, the order index n, the coefficients Cmn and
   Smn, and the uncertainties in Cmn and Smn. The (delimited) data require
   107 ASCII characters; these are followed by a pad of 13 unspecified ASCII
   characters, an ASCII carriage-return, and an ASCII line-feed."

  OBJECT                   = COLUMN
    NAME                         = "COEFFICIENT DEGREE"
    DATA_TYPE                    = ASCII_INTEGER
    START_BYTE                   = 1
    BYTES                        = 5
    FORMAT                       = "I5"
    UNIT                         = "N/A"
    DESCRIPTION                  = "The degree index m of the C and S
     coefficients in this record."
  END_OBJECT               = COLUMN

  OBJECT                   = COLUMN
    NAME                         = "COEFFICIENT ORDER"
    DATA_TYPE                    = ASCII_INTEGER
    START_BYTE                   = 7
    BYTES                        = 5
    FORMAT                       = "I5"
    UNIT                         = "N/A"
    DESCRIPTION                  = "The order index n of the C and S
     coefficients in this record."
  END_OBJECT               = COLUMN

  OBJECT                   = COLUMN
    NAME                         = "C"
    DATA_TYPE                    = ASCII_REAL
    START_BYTE                   = 13
    BYTES                        = 23
    FORMAT                       = "E23.16"
    UNIT                         = "N/A"
    DESCRIPTION                  = "The coefficient Cmn for this spherical
     harmonic model."
  END_OBJECT               = COLUMN

  OBJECT                   = COLUMN
    NAME                         = "S"
    DATA_TYPE                    = ASCII_REAL
    START_BYTE                   = 37
    BYTES                        = 23
    FORMAT                       = "E23.16"
    UNIT                         = "N/A"
    DESCRIPTION                  = "The coefficient Smn for this spherical
     harmonic model."
  END_OBJECT               = COLUMN

  OBJECT                   = COLUMN
    NAME                         = "C UNCERTAINTY"
    DATA_TYPE                    = ASCII_REAL
    START_BYTE                   = 61
    BYTES                        = 23
    FORMAT                       = "E23.16"
    UNIT                         = "N/A"
    DESCRIPTION                  = "The uncertainty in the coefficient Cmn
     for this spherical harmonic model."
  END_OBJECT               = COLUMN

  OBJECT                   = COLUMN
    NAME                         = "S UNCERTAINTY"
    DATA_TYPE                    = ASCII_REAL
    START_BYTE                   = 85
    BYTES                        = 23
    FORMAT                       = "E23.16"
    UNIT                         = "N/A"
    DESCRIPTION                  = "The uncertainty in the coefficient Smn
     for this spherical harmonic model."
  END_OBJECT               = COLUMN

END_OBJECT           = SHADR_COEFFICIENTS_TABLE

END

