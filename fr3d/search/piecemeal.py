import sys

from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from sqlalchemy.orm import aliased
from sqlalchemy.orm import session

#   Only needed for pretty printing
from sqlalchemy.dialects import mysql

from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import Column, Integer, String
from sqlalchemy import and_, or_, not_

from cred import hub


Base = declarative_base()


debug = False
debug = True
limit_val = 20


class UnitPairsInteractions(Base):
    __tablename__ = 'unit_pairs_interactions'

    unit_pairs_interactions_id = Column(Integer, primary_key=True)
    unit_id_1 = Column(String)
    unit_id_2 = Column(String)
    pdb_id = Column(String)
    f_lwbp = Column(String)
    f_stacks = Column(String)
    f_bphs = Column(String)
    f_brbs = Column(String)
    f_crossing = Column(Integer)


class UnitInfo(Base):
    __tablename__ = 'unit_info'

    unit_id = Column(String, primary_key=True)
    pdb_id = Column(String)
    model = Column(Integer)
    chain = Column(String)
    number = Column(Integer)
    unit = Column(String)
    #alt_id = Column(String)
    #ins_code = Column(String)
    #sym_op = Column(String)
    #chain_index = Column(Integer)
    #unit_type_id = Column(String)


#class SpeciesMapping(Base):
#    __tablename__ = 'species_mapping'
#
#    species_mapping_id = Column(Integer, primary_key=True)


class ExpSeqPosition(Base):
    __tablename__ = 'exp_seq_position'

    exp_seq_position_id = Column(Integer, primary_key=True)
    index = Column(Integer)
    unit = Column(String)
    exp_seq_id = Column(Integer)


class ExpSeqUnitMapping(Base):
    __tablename__= 'exp_seq_unit_mapping'

    exp_seq_unit_mapping_id = Column(Integer, primary_key=True)
    unit_id = Column(String)
    exp_seq_position_id = Column(Integer)


def conformn_terms(u1, u2, upi, families):
    return upi.f_lwbp.in_(families)


def sequence_terms(u1, u2, upi, sequences):
    terms = ((u1.unit == sequences[0][0]) & (u2.unit == sequences[0][1]))

    for seq in sequences[1:]:
        terms |= ((u1.unit == seq[0]) & (u2.unit == seq[1]))
    
    return terms


def crossing_terms(u1, u2, upi, cutoff):
    #   incomplete -- does not handle case where cutoff[0] = cutoff[1]
    return (upi.f_crossing >= cutoff[0]) & (upi.f_crossing <= cutoff[1])


def setup_db(uri):
    engine = create_engine(uri)
    return sessionmaker(bind=engine)


def pretty_print(arg):
    print(arg.statement.compile(compile_kwargs={"literal_binds": True},
                                    dialect=mysql.dialect()))


def debug_output(arg, limit_val):
    pretty_print(arg)

    for row in arg.limit(limit_val):
        print(row)


##############################################################################
#
#   JOINs (define each within function)
#
def joinNT(query, upi, u1, u2):
    return query.join(u1, upi.unit_id_1 == u1.unit_id)\
                .join(u2, upi.unit_id_2 == u2.unit_id)\
                .filter(u1.pdb_id == u2.pdb_id)\
                .filter(u1.model == u2.model)\
                .filter(u1.chain == u2.chain)


#def joinNT(query, upi, u1, u2, nt):
    #
    #   OLD VERSION (allows cross-chain interactions, which wreaks havoc
    #   with sequential distance assessments)
    #
    #if (nt == 1):
    #    return query.join(u1, upi.unit_id_1 == u1.unit_id)
    #elif (nt == 2):
    #    return query.join(u2, upi.unit_id_2 == u2.unit_id)
    #else:
    #    return query


def joinBP(query, upi, esp1, esp2, esum1, esum2):
    return query.join(esum1, upi.unit_id_1 == esum1.unit_id)\
                .join(esum2, upi.unit_id_2 == esum2.unit_id)\
                .join(esp1, esum1.exp_seq_position_id == esp1.exp_seq_position_id)\
                .join(esp2, esum2.exp_seq_position_id == esp2.exp_seq_position_id)
##############################################################################


##############################################################################
#
#   Filters
#
#   Uses the default SQL Alchemy behavior of ANDing multiple .filter calls.
#   This works as long as the implicit AND between condition types applies.
#


def parseOperation(query,upi,esp1,esp2,op,nu):
    """
    Build filter conditions for the five "Sequential Distance Constraints"
    cases.

    :op: Mathematical operator.
    :nu: Numeric comparison value.
    """
    #   N.B.  This would be simpler if we control the input so that 
    #   esp1.index is guaranteed to be less than esp2.index.  Avoids
    #   need for absolute value kludging.  Each case would reduce to 
    #   single condition (and ~22 lines of code disappear).
    if ( op == '=' ):
        line = or_(
                    ((esp1.index - esp2.index) == nu),
                    ((esp2.index - esp1.index) == nu))
    elif ( op == "<" ):
        line = or_(
                    and_(
                            ((esp1.index - esp2.index) > 0),
                            ((esp1.index - esp2.index) < nu)
                        ),
                    and_(
                            ((esp2.index - esp1.index) > 0),
                            ((esp2.index - esp1.index) < nu)
                        ))
    elif ( op == ">" ):
        line = or_(
                        ((esp1.index - esp2.index) > nu),
                        ((esp2.index - esp1.index) > nu))
    elif ( op == "<=" ):
        line = or_(
                    and_(
                            ((esp1.index - esp2.index) > 0),
                            ((esp1.index - esp2.index) <= nu)
                        ),
                    and_(
                            ((esp2.index - esp1.index) > 0),
                            ((esp2.index - esp1.index) <= nu)
                        ))
    elif ( op == ">=" ):
        line = or_(
                        ((esp1.index - esp2.index) >= nu),
                        ((esp2.index - esp1.index) >= nu))

    return line


def filterNTDistance(query, upi, esp1, esp2, distlist):
    """
    Add filter for nucleotide distance.

    May include 1 or 2 sets of criteria.  Each set contains a mathematical 
    operator and a numeric value.

    :distlist: List of criteria sets.  Each set is a list containing one 
        mathematical operator and a numeric value.
    """
    #
    #   Separate the single and double criteria cases.
    #
    #   Single cases:  five operators ( = < > <= >= )
    #   Double cases:  four operators ( < > <= >= )
    #
    #   Consider input preprocessing to put data into le/ge form.
    #   Also confirm that '=' is not an operator in a double case then.
    if (debug):
        print "list length = %s" % (len(distlist))

    op, nu = distlist[0]
    query = query.filter(parseOperation(query,upi,esp1,esp2,op,nu))
    
    if (len(distlist) == 2):
        #   AND (if present), so no extra shenanigans
        op, nu = distlist[1]
        query = query.filter(parseOperation(query,upi,esp1,esp2,op,nu))

    return query


def filterNTOrder(query, esp1, esp2, order):
    """
    Add filter for nucleotide order.

    :order: one of "<" and ">" (mutually exclusive options), given as 2 vs. 1
        via the web form
    """
    if ( order == "<" ):
        return query.filter(esp1.index > esp2.index)
    else:
        return query.filter(esp1.index < esp2.index)


def filterNTIdentity(query, upi, u1, u2, nt, nts):
    """
    Add filter for nucleotide identity to query.

    :nt: 1 or 2 (first or second nucleotide in pair).
    :nts: List of one or more nucleotides.
    """
    #   Depending upon where ambiguous nucleotides are expanded, this might
    #   become more complicated...
    if (nt == 1):
        return query.filter(u1.unit.in_(nts))
    elif (nt == 2):
        return query.filter(u2.unit.in_(nts))
    
    return query


def filterBPconf(query, upi, conflist):
    """
    Add filter for basepair conformation to query.

    :conflist: List of one or more conformations.
    """
    #   Depending upon where the category labels are expanded, this might
    #   become more complicated...
    return query.filter(upi.f_lwbp.in_(conflist))


def filterStack(query, upi, stlist):
    """
    Add filter for base stacking interactions to query.

    :stlist: List of one or more stacking patterns.
    """
    #   Depending upon where the category labels are expanded, this might
    #   become more complicated...
    return query.filter(upi.f_stacks.in_(stlist))


def filterBPHS(query, upi, bphslist):
    """
    Add filter for base-phosphate interactions to query.

    :bphslist: List of one or more base-phosphate interactions.
    """
    #   Depending upon where the category labels are expanded, this might
    #   become more complicated...
    return query.filter(upi.f_bphs.in_(bphslist))


def filterBRBS(query, upi, brbslist):
    """
    Add filter for base-ribose interactions to query.

    :brbslist: List of one or more base-ribose insteractions.
    """
    #   Depending upon where the category labels are expanded, this might
    #   become more complicated...
    return query.filter(upi.f_brbs.in_(brbslist))


def filterCrossCategory(query, upi, category):
    """
    Add filter for crossing number, based on named categories, to query.

    :category: string category name (from nested, local, distant).
    """
    #   Obviously, this does not handle "not-category" queries at present.
    #   TO DO:  process category as list (will help support NOT queries).
    if (category == 'nested'):
        cross = [0,0]
    elif (category == 'local'):
        cross = [1,2]
    elif ( category == 'distant'):
        cross = [3,0]

    return filterCrossNumeric(query, upi, cross)


def filterCrossNumeric(query, upi, cross):
    """
    Add filter for crossing number to query.

    Pre-process the input to provide the le/ge form of inequalities.
    """
    lo = cross[0]
    hi = cross[1]

    if ( hi == lo ):
        query = query.filter(upi.f_crossing == lo)
    else:
        if ( hi > lo ):
            query = query.filter(upi.f_crossing <= hi)

        if ( lo > 0 ):
            query = query.filter(upi.f_crossing >= lo)

    return query


##############################################################################


def main(db_uri, constraints):
    GENERATORS = {
        'conformn': conformn_terms,
        'sequence': sequence_terms,
        'crossing': crossing_terms
    }

    Session = setup_db(db_uri)

    sess = Session()

    u1 = aliased(UnitInfo)
    u2 = aliased(UnitInfo)

    upi = UnitPairsInteractions

    #sm = SpeciesMapping

    esp1 = aliased(ExpSeqPosition)
    esp2 = aliased(ExpSeqPosition)

    esum1 = aliased(ExpSeqUnitMapping)
    esum2 = aliased(ExpSeqUnitMapping)

    if (debug):
        for key, value in constraints.items():
            print("%s: %s") % (key, value)

    #
    #   Base query -- additional columns for debugging version
    #
    if (debug):
        query = sess.query(upi.unit_id_1, upi.unit_id_2, upi.pdb_id, 
                           upi.f_lwbp, upi.f_stacks, upi.f_bphs, 
                           upi.f_brbs, upi.f_crossing, esp1.index,
                           esp2.index)
    else:
        query = sess.query(upi.unit_id_1, upi.unit_id_2)

    #
    #   Do all JOINs here -- need the constraints between both copies
    #   of unit_info to enforce the same chain for all interactions.
    #   Also makes all columns available for debugging when needed.
    #
    query = joinNT(query, upi, u1, u2)
    query = joinBP(query, upi, esp1, esp2, esum1, esum2)

    ################################################################
    #   Hardwiring variables for testing
    #   allows for mix-match queries; with debug on, can see SQL on stdout
    #

    stack = ''

    bphs = ''
    brbs = ''

    cross = ''
    #cross = constraints['crossing']
    #cross = [2,4]
    #cross = [2,2]
    #cross = [0,3]
    #cross = [3,0]

    category = ''
    #category = 'nested'
    #category = 'local'
    #category = 'distant'

    order = ''
    #order = '<'
    #order = '>'

    nt = ''
    #nt = 1
    #nt = 2

    nts = ''
    #nts = ['A', 'C']

    distlist = ''
    #distlist = [['=', 3]]
    #distlist = [['<', 3]]
    #distlist = [['>', 20]]
    #distlist = [['<=', 3]]
    #distlist = [['>=', 20]]
    #distlist = [('>', 5), ('<=', 12)]
    ################################################################

    ################################################################
    #   Apply filters
    #   Revise to use Blake's structure -- it's cleaner.
    #
    for key, value in constraints.items():
        if key in GENERATORS:
            fn = GENERATORS[key]
            query = query.filter(fn(u1, u2, upi, value))

    if (stack):
        query = filterStack(query, upi, stack)

    if (bphs):
        query = filterBPHS(query, upi, bphs)

    if (brbs):
        query = filterBRBS(query, upi, brbs)

    if (cross):
        #   too many assumptions on the content of cross ATM
        #   should probably regularize before this step
        query = filterCrossNumeric(query, upi, cross)

    if (category):
        #   This logic won't handle NOTs correctly at present
        query = filterCrossCategory(query, upi, category)

    if (order):
        query = filterNTOrder(query, esp1, esp2, order)

    if (nt and nts):
        query = filterNTIdentity(query, upi, u1, u2, nt, nts)

    if (distlist):
        query = filterNTDistance(query, upi, esp1, esp2, distlist)

    ################################################################



    if (debug):
        debug_output(query,limit_val)


if __name__ == '__main__':
    query = {
        'conformn': ['cWW', 'tWW'],
        'sequence': ['AU', 'CG'],
        #'stack': ['s53', 's35'],
        #'stack': ['s55', 's33'],
        #'distance': ['>', 2],
        #'crossing': (1, 5)
    }
    
    main(hub, query)

quit()