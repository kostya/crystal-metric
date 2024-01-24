require "digest/crc32"

abstract class Benchmark
  abstract def run
  abstract def result
  abstract def expected

  def self.run
    summary_time = 0.0
    ok = 0
    fails = 0
    silent = ENV["SILENT"]? == "1"

    {% for a in @type.subclasses %}
      print "{{a}}: " unless silent
      bench = {{a.id}}.new
      GC.collect
      t = Time.local
      bench.run
      delta = (Time.local - t).to_f
      GC.collect
      if bench.result == bench.expected
        print "ok " unless silent
        ok += 1
      else
        print "err result=#{bench.result.inspect}, but expected=#{bench.expected.inspect} " unless silent
        fails += 1
      end
      print "in %.3fs\n" % {delta} unless silent
      summary_time += delta
      GC.collect
      sleep 0.1
      GC.collect
      sleep 0.1
    {% end %}

    unless silent
      puts "----"
      if fails > 0
        puts "ERR #{fails}"
        puts "%.4fs" % {summary_time}
        exit 1
      else
        puts "OK #{ok}"
        puts "%.4fs" % {summary_time}
      end
    else
      puts "%.4f" % {summary_time}
    end
  end
end

class Binarytrees < Benchmark
  getter result

  class TreeNode
    property left : TreeNode?
    property right : TreeNode?
    property item : Int32

    def self.create(item, depth) : TreeNode
      TreeNode.new item, depth - 1
    end

    def initialize(@item, depth = 0)
      if depth > 0
        self.left = TreeNode.new 2 * item - 1, depth - 1
        self.right = TreeNode.new 2 * item, depth - 1
      end
    end

    def check
      return item if (lft = left).nil?
      return item if (rgt = right).nil?
      lft.check - rgt.check + item
    end
  end

  def initialize(@n = 18)
    @result = 0
  end

  def run
    min_depth = 4
    max_depth = Math.max min_depth + 2, @n
    stretch_depth = max_depth + 1
    @result += TreeNode.create(0, stretch_depth).check

    long_lived_tree = TreeNode.create 0, max_depth
    min_depth.step(to: max_depth, by: 2) do |depth|
      iterations = 1 << (max_depth - depth + min_depth)
      1.upto(iterations) do |i|
        @result += TreeNode.create(i, depth).check
        @result += TreeNode.create(-i, depth).check
      end
    end
  end

  def expected
    -699041
  end
end

class Brainfuck < Benchmark
  getter result

  class Tape
    def initialize
      @tape = [0]
      @pos = 0
    end

    @[AlwaysInline]
    def get
      @tape[@pos]
    end

    @[AlwaysInline]
    def inc
      @tape[@pos] += 1
    end

    @[AlwaysInline]
    def dec
      @tape[@pos] -= 1
    end

    @[AlwaysInline]
    def advance
      @pos += 1
      @tape << 0 if @tape.size <= @pos
    end

    @[AlwaysInline]
    def devance
      @pos -= 1 if @pos > 0
    end
  end

  class Program
    def initialize(text)
      @chars = [] of Char
      @bracket_map = {} of Int32 => Int32
      leftstack = [] of Int32
      pc = 0
      text.each_char do |char|
        if "[]<>+-,.".includes?(char)
          @chars << char
          if char == '['
            leftstack << pc
          elsif char == ']' && !leftstack.empty?
            left = leftstack.pop
            right = pc
            @bracket_map[left] = right
            @bracket_map[right] = left
          end
          pc += 1
        end
      end
    end

    def run
      result = 0
      tape = Tape.new
      pc = 0
      while pc < @chars.size
        case @chars[pc]
        when '+'; tape.inc
        when '-'; tape.dec
        when '>'; tape.advance
        when '<'; tape.devance
        when '['; pc = @bracket_map[pc] if tape.get == 0
        when ']'; pc = @bracket_map[pc] if tape.get != 0
        when '.'; result += tape.get.chr.ord
        else
        end
        pc += 1
      end
      result
    end
  end

  def initialize
    @text = "Benchmark brainf*ck program
>++[<+++++++++++++>-]<[[>+>+<<-]>[<+>-]++++++++
[>++++++++<-]>.[-]<<>++++++++++[>++++++++++[>++
++++++++[>++++++++++[>++++++++++[>++++++++++[>+
+++++++++[-]<-]<-]<-]<-]<-]<-]<-]++++++++++."
    @result = 0
  end

  def run
    @result = Program.new(@text).run
  end

  def expected
    2025
  end
end

class Brainfuck2 < Benchmark
  getter result

  module Op
    record Inc, val : Int32
    record Move, val : Int32
    record Print
    alias T = Inc | Move | Print | Array(Op::T)
  end

  class Tape
    def initialize
      @tape = [0_u8]
      @pos = 0
    end

    @[AlwaysInline]
    def get
      @tape[@pos]
    end

    @[AlwaysInline]
    def inc(x)
      @tape[@pos] += x
    end

    @[AlwaysInline]
    def move(x)
      @pos += x
      while @pos >= @tape.size
        @tape << 0
      end
    end
  end

  class Program
    property result
    @ops : Array(Op::T)

    def initialize(code : String)
      @ops = parse(code.each_char)
      @result = 0
    end

    def run
      _run @ops, Tape.new
    end

    private def _run(program, tape)
      program.each do |op|
        case op
        when Op::Inc
          tape.inc(op.val)
        when Op::Move
          tape.move(op.val)
        when Array(Op::T)
          while tape.get != 0
            _run(op, tape)
          end
        when Op::Print
          @result += tape.get.chr.ord
        else
        end
      end
    end

    private def parse(iterator)
      res = [] of Op::T
      iterator.each do |c|
        op = case c
             when '+'; Op::Inc.new(1)
             when '-'; Op::Inc.new(-1)
             when '>'; Op::Move.new(1)
             when '<'; Op::Move.new(-1)
             when '.'; Op::Print.new
             when '['; parse(iterator)
             when ']'; break
             else
             end
        res << op if op
      end
      res
    end
  end

  def initialize
    @text = "Benchmark brainf*ck program
>++[<+++++++++++++>-]<[[>+>+<<-]>[<+>-]++++++++
[>++++++++<-]>.[-]<<>++++++++++[>++++++++++[>++
++++++++[>++++++++++[>++++++++++[>++++++++++[>+
+++++++++[-]<-]<-]<-]<-]<-]<-]<-]++++++++++."
    @result = 0
  end

  def run
    program = Program.new(@text)
    program.run
    @result = program.result
  end

  def expected
    2025
  end
end

class Fannkuchredux < Benchmark
  getter result

  def fannkuchredux(n)
    perm1 = StaticArray(Int32, 32).new { |i| i }
    perm = StaticArray(Int32, 32).new(0)
    count = StaticArray(Int32, 32).new(0)
    maxFlipsCount = permCount = checksum = 0
    r = n

    while true
      while r > 1
        count[r - 1] = r
        r -= 1
      end

      n.times { |i| perm[i] = perm1[i] } # dup
      flipsCount = 0

      while !((k = perm[0]) == 0)
        k2 = (k + 1) >> 1
        (0...k2).each do |i|
          j = k - i
          perm[i], perm[j] = perm[j], perm[i] # swap
        end
        flipsCount += 1
      end

      maxFlipsCount = flipsCount if flipsCount > maxFlipsCount
      checksum += (permCount % 2 == 0) ? flipsCount : -flipsCount

      while true
        return {checksum, maxFlipsCount} if r == n

        perm0 = perm1[0]
        (0...r).each do |i|
          j = i + 1
          perm1[i], perm1[j] = perm1[j], perm1[i] # swap
        end

        perm1[r] = perm0
        cntr = count[r] -= 1
        break if cntr > 0
        r += 1
      end
      permCount += 1
    end
  end

  def initialize(@n = 11)
    @result = {0, 0}
  end

  def run
    @result = fannkuchredux(@n)
  end

  def expected
    {556355, 51}
  end
end

class Fasta < Benchmark
  IM = 139968
  IA =   3877
  IC =  29573

  class Last
    def self.last
      @@last ||= 42
    end

    def self.last=(x)
      @@last = x
    end
  end

  def gen_random(max)
    Last.last = (Last.last * IA + IC) % IM
    max * Last.last / IM.to_f64
  end

  def make_cumulative(genelist)
    cp = 0.0_f64
    genelist.size.times do |i|
      c, p = genelist[i]
      cp += p
      genelist[i] = {c, cp}
    end
  end

  def select_random(genelist)
    r = gen_random(1)
    return genelist[0][0] if r < genelist[0][1]

    lo = 0
    hi = genelist.size - 1

    while hi > lo + 1
      i = (hi + lo) // 2
      if r < genelist[i][1]
        hi = i
      else
        lo = i
      end
    end
    genelist[hi][0]
  end

  LINE_LENGTH = 60

  def make_random_fasta(id, desc, genelist, n)
    todo = n
    @result << ">#{id} #{desc}"
    @result << "\n"

    while todo > 0
      m = (todo < LINE_LENGTH) ? todo : LINE_LENGTH
      pick = String.new(m) do |buffer|
        m.times { |i| buffer[i] = select_random(genelist).ord.to_u8 }
        {m, m}
      end
      @result << pick
      @result << "\n"
      todo -= LINE_LENGTH
    end
  end

  def make_repeat_fasta(id, desc, s, n)
    todo = n
    k = 0
    kn = s.size

    @result << ">#{id} #{desc}"
    @result << "\n"
    while todo > 0
      m = (todo < LINE_LENGTH) ? todo : LINE_LENGTH

      while m >= kn - k
        @result << s[k..-1]
        m -= kn - k
        k = 0
      end

      @result << s[k...k + m]
      @result << "\n"
      k += m

      todo -= LINE_LENGTH
    end
  end

  IUB = [
    {'a', 0.27},
    {'c', 0.12},
    {'g', 0.12},
    {'t', 0.27},

    {'B', 0.02},
    {'D', 0.02},
    {'H', 0.02},
    {'K', 0.02},
    {'M', 0.02},
    {'N', 0.02},
    {'R', 0.02},
    {'S', 0.02},
    {'V', 0.02},
    {'W', 0.02},
    {'Y', 0.02},
  ]

  HOMO = [{'a', 0.3029549426680}, {'c', 0.1979883004921}, {'g', 0.1975473066391}, {'t', 0.3015094502008}]
  ALU  = "GGCCGGGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGGCGGGCGGATCACCTGAGGTCAGGAGTTCGAGACCAGCCTGGCCAACATGGTGAAACCCCGTCTCTACTAAAAATACAAAAATTAGCCGGGCGTGGTGGCGCGCGCCTGTAATCCCAGCTACTCGGGAGGCTGAGGCAGGAGAATCGCTTGAACCCGGGAGGCGGAGGTTGCAGTGAGCCGAGATCGCGCCACTGCACTCCAGCCTGGGCGACAGAGCGAGACTCCGTCTCAAAAA"

  def initialize(@n = 7000000)
    make_cumulative(IUB)
    make_cumulative(HOMO)
    @result = IO::Memory.new
  end

  def run
    make_repeat_fasta("ONE", "Homo sapiens alu", ALU, @n * 2)
    make_random_fasta("TWO", "IUB ambiguity codes", IUB, @n * 3)
    make_random_fasta("THREE", "Homo sapiens frequency", HOMO, @n * 5)
  end

  def result
    Digest::CRC32.checksum(@result.to_s)
  end

  def expected
    223599119
  end
end

class Knuckeotide < Benchmark
  @seq : String

  def frecuency(seq, length)
    n = seq.size - length + 1
    table = Hash(String, Int32).new { 0 }
    (0...n).each do |f|
      table[seq.byte_slice(f, length)] += 1
    end
    {n, table}
  end

  def sort_by_freq(seq, length)
    n, table = frecuency(seq, length)
    table.to_a.sort { |a, b| b[1] <=> a[1] }.each do |v|
      @result << "%s %.3f\n" % {v[0].upcase, ((v[1] * 100).to_f / n)}
    end
    @result << "\n"
  end

  def find_seq(seq, s)
    n, table = frecuency(seq, s.size)
    @result << "#{table[s].to_s}\t#{s.upcase}\n"
  end

  def run(seq)
  end

  def initialize
    @result = IO::Memory.new

    f = Fasta.new(600000)
    f.run
    res = f.@result.to_s

    seq = IO::Memory.new
    three = false

    res.each_line do |line|
      if line.starts_with?(">THREE")
        three = true
        next
      end
      seq << line.chomp if three
    end
    @seq = seq.to_s
  end

  def run
    (1..2).each { |i| sort_by_freq(@seq, i) }
    %w(ggt ggta ggtatt ggtattttaatt ggtattttaatttatagt).each { |s| find_seq(@seq, s) }
  end

  def result
    Digest::CRC32.checksum(@result.to_s)
  end

  def expected
    159474310
  end
end

class Mandelbrot < Benchmark
  ITER  =  50
  LIMIT = 2.0

  def initialize(@n = 3800)
    @result = IO::Memory.new
  end

  def run
    w = h = @n
    @result << "P4\n#{w} #{h}\n"

    bit_num = 0
    byte_acc = 0_u8

    h.times do |y|
      w.times do |x|
        zr = zi = tr = ti = 0.0
        cr = (2.0 * x / w - 1.5)
        ci = (2.0 * y / h - 1.0)

        i = 0
        while (i < ITER) && (tr + ti <= LIMIT * LIMIT)
          zi = 2.0 * zr * zi + ci
          zr = tr - ti + cr
          tr = zr * zr
          ti = zi * zi
          i += 1
        end

        byte_acc <<= 1
        byte_acc |= 0x01 if tr + ti <= LIMIT * LIMIT
        bit_num += 1

        if bit_num == 8
          @result.write_byte byte_acc
          byte_acc = 0_u8
          bit_num = 0
        elsif x == w - 1
          byte_acc <<= 8 - w % 8
          @result.write_byte byte_acc
          byte_acc = 0_u8
          bit_num = 0
        end
      end
    end
  end

  def result
    Digest::CRC32.checksum(@result.to_s)
  end

  def expected
    2862550248
  end
end

class Matmul < Benchmark
  def matmul(a, b)
    m = a.size
    n = a[0].size
    p = b[0].size
    # transpose
    b2 = Array.new(n) { Array.new(p, 0.0) }
    (0...n).each do |i|
      (0...p).each do |j|
        b2[j][i] = b[i][j]
      end
    end
    # multiplication
    c = Array.new(m) { Array.new(p, 0.0) }
    c.each_with_index do |ci, i|
      ai = a[i]
      b2.each_with_index do |b2j, j|
        s = 0.0
        b2j.each_with_index do |b2jv, k|
          s += ai[k] * b2jv
        end
        ci[j] = s
      end
    end
    c
  end

  def matgen(n)
    tmp = 1.0 / n / n
    a = Array.new(n) { Array.new(n, 0.0) }
    (0...n).each do |i|
      (0...n).each do |j|
        a[i][j] = tmp * (i - j) * (i + j)
      end
    end
    a
  end

  getter result

  def initialize(@n = 980)
    @result = 0.0
  end

  def run
    a = matgen(@n)
    b = matgen(@n)
    c = matmul(a, b)
    @result = c[@n >> 1][@n >> 1]
  end

  def expected
    -93.66692176867205
  end
end

class Nbody < Benchmark
  # Copied with little modifications from: http://benchmarksgame.alioth.debian.org/u32/program.php?test=nbody&lang=yarv&id=2
  SOLAR_MASS    = 4 * Math::PI**2
  DAYS_PER_YEAR = 365.24

  class Planet
    def_clone

    property x : Float64
    property y : Float64
    property z : Float64
    property vx : Float64
    property vy : Float64
    property vz : Float64
    property mass : Float64

    def initialize(@x, @y, @z, vx, vy, vz, mass)
      @vx, @vy, @vz = vx * DAYS_PER_YEAR, vy * DAYS_PER_YEAR, vz * DAYS_PER_YEAR
      @mass = mass * SOLAR_MASS
    end

    def move_from_i(bodies, nbodies, dt, i)
      while i < nbodies
        b2 = bodies[i]
        dx = @x - b2.x
        dy = @y - b2.y
        dz = @z - b2.z

        distance = Math.sqrt(dx * dx + dy * dy + dz * dz)
        mag = dt / (distance * distance * distance)
        b_mass_mag, b2_mass_mag = @mass * mag, b2.mass * mag

        @vx -= dx * b2_mass_mag
        @vy -= dy * b2_mass_mag
        @vz -= dz * b2_mass_mag
        b2.vx += dx * b_mass_mag
        b2.vy += dy * b_mass_mag
        b2.vz += dz * b_mass_mag
        i += 1
      end

      @x += dt * @vx
      @y += dt * @vy
      @z += dt * @vz
    end
  end

  def energy(bodies)
    e = 0.0
    nbodies = bodies.size

    0.upto(nbodies - 1) do |i|
      b = bodies[i]
      e += 0.5 * b.mass * (b.vx * b.vx + b.vy * b.vy + b.vz * b.vz)
      (i + 1).upto(nbodies - 1) do |j|
        b2 = bodies[j]
        dx = b.x - b2.x
        dy = b.y - b2.y
        dz = b.z - b2.z
        distance = Math.sqrt(dx * dx + dy * dy + dz * dz)
        e -= (b.mass * b2.mass) / distance
      end
    end
    e
  end

  def offset_momentum(bodies)
    px, py, pz = 0.0, 0.0, 0.0

    bodies.each do |b|
      m = b.mass
      px += b.vx * m
      py += b.vy * m
      pz += b.vz * m
    end

    b = bodies[0]
    b.vx = -px / SOLAR_MASS
    b.vy = -py / SOLAR_MASS
    b.vz = -pz / SOLAR_MASS
  end

  BODIES = [
    # sun
    Planet.new(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0),

    # jupiter
    Planet.new(
      4.84143144246472090e+00,
      -1.16032004402742839e+00,
      -1.03622044471123109e-01,
      1.66007664274403694e-03,
      7.69901118419740425e-03,
      -6.90460016972063023e-05,
      9.54791938424326609e-04),

    # saturn
    Planet.new(
      8.34336671824457987e+00,
      4.12479856412430479e+00,
      -4.03523417114321381e-01,
      -2.76742510726862411e-03,
      4.99852801234917238e-03,
      2.30417297573763929e-05,
      2.85885980666130812e-04),

    # uranus
    Planet.new(
      1.28943695621391310e+01,
      -1.51111514016986312e+01,
      -2.23307578892655734e-01,
      2.96460137564761618e-03,
      2.37847173959480950e-03,
      -2.96589568540237556e-05,
      4.36624404335156298e-05),

    # neptune
    Planet.new(
      1.53796971148509165e+01,
      -2.59193146099879641e+01,
      1.79258772950371181e-01,
      2.68067772490389322e-03,
      1.62824170038242295e-03,
      -9.51592254519715870e-05,
      5.15138902046611451e-05),
  ]

  def initialize(@n = 15000000)
    @result = {0.0, 0.0}
    @body = BODIES
  end

  def run
    offset_momentum(@body)

    v1 = energy(@body)

    nbodies = @body.size
    dt = 0.01

    @n.times do
      i = 0
      while i < nbodies
        b = @body[i]
        b.move_from_i(@body, nbodies, dt, i + 1)
        i += 1
      end
    end

    v2 = energy(@body)
    @result = {v1, v2}
  end

  getter result

  def expected
    {-0.16907516382852447, -0.16905678767532126}
  end
end

require "big"

class Pidigits < Benchmark
  def initialize(@nn = 4500)
    @result = IO::Memory.new
  end

  def run
    start = Time.local

    i = 0
    k = 0
    ns = 0.to_big_i
    a = 0.to_big_i
    t = 0
    u = 0.to_big_i
    k1 = 1
    n = 1.to_big_i
    d = 1.to_big_i

    while true
      k += 1
      t = n << 1
      n *= k
      k1 += 2
      a = (a + t) * k1
      d *= k1
      if a >= n
        t, u = (n * 3 + a).divmod(d)
        u += n
        if d > u
          ns = ns * 10 + t
          i += 1
          if i % 10 == 0
            @result << "%010d\t:%d\n" % {ns.to_u64, i}
            ns = 0
          end
          break if i >= @nn

          a = (a - (d * t)) * 10
          n *= 10
        end
      end
    end
  end

  def result
    Digest::CRC32.checksum(@result.to_s)
  end

  def expected
    2259604824
  end
end

class RegexDna < Benchmark
  @seq : String
  @ilen : Int32
  @clen : Int32

  def initialize
    @result = IO::Memory.new

    f = Fasta.new(6000000)
    f.run
    res = f.@result.to_s

    seq = IO::Memory.new

    @ilen = 0
    res.each_line do |line|
      @ilen += line.bytesize + 1
      seq << line.chomp unless line.starts_with? '>'
    end

    @seq = seq.to_s
    @clen = seq.bytesize
  end

  def run
    [
      /agggtaaa|tttaccct/,
      /[cgt]gggtaaa|tttaccc[acg]/,
      /a[act]ggtaaa|tttacc[agt]t/,
      /ag[act]gtaaa|tttac[agt]ct/,
      /agg[act]taaa|ttta[agt]cct/,
      /aggg[acg]aaa|ttt[cgt]ccct/,
      /agggt[cgt]aa|tt[acg]accct/,
      /agggta[cgt]a|t[acg]taccct/,
      /agggtaa[cgt]|[acg]ttaccct/,
    ].each { |f| @result << "#{f.source} #{@seq.scan(f).size}\n" }

    hash = {
      "B" => "(c|g|t)",
      "D" => "(a|g|t)",
      "H" => "(a|c|t)",
      "K" => "(g|t)",
      "M" => "(a|c)",
      "N" => "(a|c|g|t)",
      "R" => "(a|g)",
      "S" => "(c|t)",
      "V" => "(a|c|g)",
      "W" => "(a|t)",
      "Y" => "(c|t)",
    }

    t = Time.local
    @seq = @seq.gsub(/B|D|H|K|M|N|R|S|V|W|Y/, hash)

    @result << "\n"
    @result << "#{@ilen}\n"
    @result << "#{@clen}\n"
    @result << "#{@seq.size}\n"
  end

  def result
    Digest::CRC32.checksum(@result.to_s)
  end

  def expected
    3506965167
  end
end

class Revcomp < Benchmark
  @input : String

  def revcomp(seq)
    seq = seq.reverse.tr("wsatugcyrkmbdhvnATUGCYRKMBDHVN", "WSTAACGRYMKVHDBNTAACGRYMKVHDBN")
    stringlen = seq.size - 1
    0.step(to: stringlen, by: 60) { |x| @result << seq[x...x + 60]; @result << "\n" }
  end

  def initialize
    @result = IO::Memory.new

    f = Fasta.new(13000000)
    f.run
    @input = f.@result.to_s
  end

  def run
    seq = IO::Memory.new

    @input.each_line do |line|
      if line.starts_with? '>'
        if !seq.empty?
          revcomp(seq.to_s)
          seq.clear
        end
        @result << line
        @result << "\n"
      else
        seq << line.chomp
      end
    end
    revcomp(seq.to_s)
  end

  def result
    Digest::CRC32.checksum(@result.to_s)
  end

  def expected
    3936352977
  end
end

class Spectralnorm < Benchmark
  def eval_A(i, j)
    1.0_f64 / ((i + j) * (i + j + 1.0) / 2.0 + i + 1.0)
  end

  def eval_A_times_u(u)
    (0...u.size).map do |i|
      v = 0.0_f64
      u.each_with_index do |uu, j|
        v += eval_A(i, j) * uu
      end
      v
    end
  end

  def eval_At_times_u(u)
    (0...u.size).map do |i|
      v = 0.0_f64
      u.each_with_index do |uu, j|
        v += eval_A(j, i) * uu
      end
      v
    end
  end

  def eval_AtA_times_u(u)
    eval_At_times_u(eval_A_times_u(u))
  end

  getter result

  def initialize(@n = 4500)
    @result = 0.0
  end

  def run
    u = Array.new(@n, 1.0_f64)
    v = Array.new(@n, 1.0_f64)
    10.times do
      v = eval_AtA_times_u(u)
      u = eval_AtA_times_u(v)
    end
    vBv = vv = 0.0_f64
    (0...@n).each do |i|
      vBv += u[i] * v[i]
      vv += v[i] * v[i]
    end
    @result = Math.sqrt(vBv / vv)
  end

  def expected
    1.2742241527688491
  end
end

class Threadring < Benchmark
  THREAD_COUNT = 503

  class Receiver
    def initialize(@name : Int32, @res : Channel(Int32))
      @mailbox = Channel(Int32).new
    end

    def next=(n : Receiver)
      @next = n
    end

    def put(msg)
      @mailbox.send(msg)
    end

    def messageloop
      while true
        msg = @mailbox.receive
        if msg == 0
          @res.send(@name)
        elsif nxt = @next
          nxt.put(msg - 1)
        end
      end
    end
  end

  @receivers : Array(Receiver)

  def initialize(@n = 15000000)
    @res = Channel(Int32).new
    @receivers = Array.new(THREAD_COUNT) { |i| Receiver.new(i + 1, @res) }
    (0...THREAD_COUNT - 1).each { |i| @receivers[i].next = @receivers[i + 1] }
    @receivers[THREAD_COUNT - 1].next = @receivers[0]
    @result = 0
  end

  def run
    THREAD_COUNT.times do |i|
      spawn { @receivers[i].messageloop }
    end
    @receivers[0].put(@n)
    @result = @res.receive
  end

  getter result : Int32

  def expected
    38
  end
end

require "base64"

class Base64Encode < Benchmark
  TRIES = 8192

  @str : String
  @str2 : String
  @str3 : String

  def initialize(@n = 300_000)
    @str = "a" * @n
    @str2 = Base64.strict_encode(@str)
    @str3 = Base64.decode_string(@str2)
    @result = ""
  end

  def run
    s_encoded = 0_i64

    TRIES.times do |i|
      s_encoded += Base64.strict_encode(@str).bytesize
    end

    @result += "encode #{@str[0..3]}... to #{@str2[0..3]}...: #{s_encoded}\n"
  end

  getter result

  def expected
    "encode aaaa... to YWFh...: 3276800000\n"
  end
end

class Base64Decode < Benchmark
  TRIES = 8192

  @str : String
  @str2 : String
  @str3 : String

  def initialize(@n = 230_000)
    @str = "a" * @n
    @str2 = Base64.strict_encode(@str)
    @str3 = Base64.decode_string(@str2)
    @result = ""
  end

  def run
    s_decoded = 0_i64

    TRIES.times do |i|
      s_decoded += Base64.decode_string(@str2).bytesize
    end

    @result += "decode #{@str2[0..3]}... to #{@str3[0..3]}...: #{s_decoded}\n"
  end

  getter result

  def expected
    "decode YWFh... to aaaa...: 1884160000\n"
  end
end

class Primes < Benchmark
  PREFIX = 32_338

  class Node
    property :children, :terminal

    def initialize
      @children = Hash(Char, Node).new
      @terminal = false
    end
  end

  class Sieve
    def initialize(limit : Int32)
      @limit = limit
      @prime = Array(Bool).new(limit + 1, false)
    end

    def to_list
      result = [2, 3]
      (5..@limit).each do |p|
        result.push(p) if @prime[p]
      end
      result
    end

    def omit_squares
      r = 5
      while r * r < @limit
        if @prime[r]
          i = r * r
          while i < @limit
            @prime[i] = false
            i += r * r
          end
        end
        r += 1
      end

      self
    end

    def step1(x, y)
      n = (4 * x * x) + (y * y)
      @prime[n] = !@prime[n] if n <= @limit && (n % 12 == 1 || n % 12 == 5)
    end

    def step2(x, y)
      n = (3 * x * x) + (y * y)
      @prime[n] = !@prime[n] if n <= @limit && n % 12 == 7
    end

    def step3(x, y)
      n = (3 * x * x) - (y * y)
      @prime[n] = !@prime[n] if x > y && n <= @limit && n % 12 == 11
    end

    def loop_y(x)
      y = 1
      while y * y < @limit
        step1(x, y)
        step2(x, y)
        step3(x, y)
        y += 1
      end
    end

    def loop_x
      x = 1
      while x * x < @limit
        loop_y(x)
        x += 1
      end
    end

    def calc
      loop_x
      omit_squares
    end
  end

  def generate_trie(l)
    root = Node.new
    l.each do |el|
      head = root
      el.to_s.each_char do |ch|
        head.children[ch] = Node.new unless head.children[ch]?
        head = head.children[ch]
      end
      head.terminal = true
    end
    root
  end

  def find(upper_bound, prefix)
    primes = Sieve.new(upper_bound).calc
    str_prefix = prefix.to_s
    head = generate_trie(primes.to_list)
    str_prefix.each_char do |ch|
      head = head.children[ch]
      return nil if head.nil?
    end

    queue = [{head, str_prefix}]
    result = Array(Int32).new
    until queue.empty?
      top, prefix = queue.pop
      result.push(prefix.to_i) if top.terminal
      top.children.each do |ch, v|
        queue.insert(0, {v, prefix + ch})
      end
    end
    result
  end

  @results : Array(Int32)?

  def initialize(@n = 40_000_000)
    @results = Array(Int32).new
  end

  def run
    @results = find(@n, PREFIX)
  end

  def result
    @results
  end

  def expected
    [323381, 323383, 3233803, 3233809, 3233851, 3233863, 3233873, 3233887, 3233897, 32338123, 32338129, 32338147, 32338169, 32338351, 32338367, 32338387, 32338039, 32338099, 32338013, 32338027, 32338043, 32338049, 32338067, 32338517, 32338507, 32338541, 32338577, 32338633, 32338637, 32338609, 32338619, 32338643, 32338661, 32338687, 32338697, 32338699, 32338739, 32338751, 32338763, 32338871, 32338819, 32338829, 32338837, 32338843, 32338849, 32338853, 32338861, 32338907, 32338909, 32338937, 32338949, 32338993, 32338997, 32338223, 32338231, 32338249, 32338291, 32338433, 32338499]
  end
end

require "json"

class JsonGenerate < Benchmark
  struct Coordinate
    include JSON::Serializable

    def initialize(@x : Float64, @y : Float64, @z : Float64, @name : String, @opts : Hash(String, Tuple(Int32, Bool)))
    end
  end

  def initialize(@n = 1_800_000)
    @result = IO::Memory.new
    @data = Array(Coordinate).new
    r = Random.new(0, 0_u64)
    letters = ('a'..'z').to_a
    @n.times do
      @data << Coordinate.new(
        r.next_float,
        r.next_float,
        r.next_float,
        "#{letters.sample(6, r).join} #{r.rand(10_000)}",
        {"1" => {1, true}},
      )
    end
  end

  def run
    {"coordinates": @data,
     "info":        "some info"}.to_json(@result)
    true
  end

  def result
    Digest::CRC32.checksum(@result.to_s)
  end

  def expected
    2191513974
  end
end

class JsonParsePure < Benchmark
  struct Coordinate
    property x, y, z

    def initialize(@x : Float64, @y : Float64, @z : Float64)
    end
  end

  def calc(text)
    jobj = JSON.parse(text)
    coordinates = jobj["coordinates"].as_a
    len = coordinates.size
    x = y = z = 0

    coordinates.each do |coord|
      x += coord["x"].as_f
      y += coord["y"].as_f
      z += coord["z"].as_f
    end

    Coordinate.new(x / len, y / len, z / len)
  end

  getter result
  @text : String

  def initialize
    j = JsonGenerate.new(600_000)
    j.run
    @text = j.@result.to_s
    @result = Coordinate.new(0, 0, 0)
  end

  def run
    @result = calc(@text)
  end

  def expected
    Coordinate.new(0.4999535894712113, 0.5000557786698896, 0.5007583687182091)
  end
end

class JsonParseSerializable < Benchmark
  struct Coordinate
    include JSON::Serializable

    property x : Float64
    property y : Float64
    property z : Float64

    def initialize(@x, @y, @z)
    end
  end

  class Coordinates
    include JSON::Serializable

    property coordinates : Array(Coordinate)
  end

  def calc(text)
    coordinates = Coordinates.from_json(text).coordinates
    len = coordinates.size
    x = y = z = 0

    coordinates.each do |e|
      x += e.x
      y += e.y
      z += e.z
    end

    Coordinate.new(x / len, y / len, z / len)
  end

  getter result
  @text : String

  def initialize
    j = JsonGenerate.new(800_000)
    j.run
    @text = j.@result.to_s
    @result = Coordinate.new(0.0, 0.0, 0.0)
  end

  def run
    @result = calc(@text)
  end

  def expected
    Coordinate.new(0.49986260081788564, 0.5001137023627188, 0.5006048890910882)
  end
end

class JsonParsePull < Benchmark
  struct Coordinate
    property x, y, z

    def initialize(@x : Float64, @y : Float64, @z : Float64)
    end
  end

  def calc(text)
    len = 0
    x = y = z = 0

    pull = JSON::PullParser.new(text)
    pull.on_key!("coordinates") do
      pull.read_array do
        len += 1
        pull.read_object do |key|
          case key
          when "x" then x += pull.read_float
          when "y" then y += pull.read_float
          when "z" then z += pull.read_float
          else          pull.skip
          end
        end
      end
    end

    Coordinate.new(x / len, y / len, z / len)
  end

  getter result
  @text : String

  def initialize
    j = JsonGenerate.new(800_000)
    j.run
    @text = j.@result.to_s
    @result = Coordinate.new(0, 0, 0)
  end

  def run
    @result = calc(@text)
  end

  def expected
    Coordinate.new(0.49986260081788564, 0.5001137023627188, 0.5006048890910882)
  end
end

Benchmark.run
