from django.core.management.base import BaseCommand

from lab.models import Sample, Metadata


class Command(BaseCommand):
    help = "Print ent_date from Metadata for given sample pseudonymized IDs"

    def handle(self, *args, **options):
        # 👇 put your sample IDs here
        sample_ids = ['m2-SrqwxQ', 'm2-Kg3Fdm', 'm2-mzMRtV', 'm2-wfEdma', 'm2-NDEjHn', 'm2-4L3pXZ', 'm2-F2HifQ',
                      'm2-9eyPAv', 'm2-NmHGDa', 'm2-HXHQAT', 'm2-gnDemF', 'm2-u89dBs', 'm2-3oZJ5x', 'm2-43Vmvj',
                      'm2-zwqknf', 'm2-w9v3Qv', 'm2-EJ7fEz', 'm2-Y4rZ7q', 'm2-XHdiUs', 'm2-GRjNzo', 'm2-Dq37e8',
                      'm2-fhoncL', 'm2-h2ZJd5', 'm2-SpXNmN', 'm2-K2udH5', 'm2-4t7Ncr', 'm2-hzNQRG', 'm2-uTpHe5',
                      'm2-BnbHLr', 'm2-LosGHn', 'm2-LCBdav', 'm2-sr4bvM', 'm2-iDqyV9', 'm2-XMdsys', 'm2-Y4wzR3',
                      'm2-tDscqy', 'm2-mmYEfJ', 'm2-Rtjhba', 'm2-VT3bcF', 'm2-fafZrX', 'm2-DgRuyg', 'm2-YBrMoY',
                      'm2-hozdkG', 'm2-kNU7js', 'm2-YowhiU', 'm2-gsuEGK', 'm2-phtHeq', 'm2-AcEG8S', 'm2-HtiJwT',
                      'm2-C983Zw', 'm2-JFKsEJ', 'm2-jHbrsp', 'm2-PGoRgN', 'm2-YYRzxL', 'm2-4Q8Stt', 'm2-BbET5x',
                      'm2-NESNdH', 'm2-LYceeY', 'm2-os283P', 'm2-QqAEGg', 'm2-gF3E5u', 'm2-sAyJww', 'm2-b99ooj',
                      'm2-FPVUGC', 'm2-2hKeW2', 'm2-j3ZXpL', 'm2-TPsF65', 'm2-uSTYZ2', 'm2-f8rfBe', 'm2-FXFt6J',
                      'm2-U2wv8C', 'm2-6z2wf2', 'm2-5cJWyf', 'm2-sFFUfJ', 'm2-3p7s5a', 'm2-ztFGoy', 'm2-qeHUBn',
                      'm2-bpDXvz', 'm2-wVBQbT', 'm2-33tHNX', 'm2-qWgDtq', 'm2-8ochGa', 'm2-4SnyCD', 'm2-kwbwf8',
                      'm2-UqhsDL', 'm2-RyXauM', 'm2-ydSpNQ', 'm2-SKJgwA', 'm2-DwiwwU', 'm2-xYnaUy', 'm2-xrRPFj',
                      'm2-sxkK9R', 'm2-DMLKnU', 'm2-2b8YpQ', 'm2-gyUNUP', 'm2-PZQaiq', 'm2-hPpCR3', 'm2-WrGEQp',
                      'm2-yhE88h', 'm2-gaQkyT', 'm2-TPgies', 'm2-3qDrEd', 'm2-XEuJCy', 'm2-ZAyjrq', 'm2-mc3GSV',
                      'm2-CLwScw', 'm2-4UF9K6', 'm2-jtjCjd', 'm2-4Rt9sw', 'm2-3URpmc', 'm2-3MgVni', 'm2-BdtP8F',
                      'm2-32WNFL', 'm2-3edUBd', 'm2-375EUk', 'm2-xwkVbK', 'm2-itJLqU', 'm2-xxAxC9', 'm2-g33df2',
                      'm2-opfXjQ', 'm2-Sm9GYA', 'm2-eAr6Lo', 'm2-xHGR9u', 'm2-4rvG7f', 'm2-SQZYHn', 'm2-oLXZ9a',
                      'm2-9HhvMh', 'm2-FHPnrx', 'm2-URVR7h', 'm2-ss4Fh8', 'm2-VA9Yi5', 'm2-qtZxkz', 'm2-PGaVE8',
                      'm2-39kc9U', 'm2-q2oXnd', 'm2-x6NRYc', 'm2-MHWQwo', 'm2-kBgpZ4', 'm2-voFSf6', 'm2-Qb4JVo',
                      'm2-V5iihT', 'm2-m6243p', 'm2-3VC6ft', 'm2-gkcF5C', 'm2-3H4uLZ', 'm2-juVwzj', 'm2-4H6FKQ',
                      'm2-DxMevL', 'm2-5xXjjp', 'm2-PVnipc', 'm2-CcozBt', 'm2-HUBpy4', 'm2-KvvZnz', 'm2-WpJqxj',
                      'm2-Yg94Uo', 'm2-STGkhE', 'm2-LTGmFn', 'm2-afGSAT', 'm2-pzbDPq', 'm2-zhrVRY', 'm2-Qf8cLC',
                      'm2-GbMzJy', 'm2-yebiMQ', 'm2-GmB4Bf', 'm2-FAvEvi', 'm2-k9vgiw', 'm2-HkLtpW', 'm2-cyCNuA',
                      'm2-F34WtX', 'm2-q4bsMd', 'm2-KaJ5Pe', 'm2-AqTpPL', 'm2-P9mWJr', 'm2-gJTbaX', 'm2-cWAMxe',
                      'm2-Ki8Aa3', 'm2-3eGjGj', 'm2-BdNxPR', 'm2-VAEUJ5', 'm2-3Vw5ij', 'm2-T5VXtD', 'm2-4SRYNi',
                      'm2-xPAa2d', 'm2-v3KDUF', 'm2-tUrV5z', 'm2-QRRsiL', 'm2-zGrf4n', 'm2-UZk6rp', 'm2-mVqwfv',
                      'm2-48cAtM', 'm2-5at3hA', 'm2-3G3Doj', 'm2-fngLcf', 'm2-wM9MFc', 'm2-VduGeq', 'm2-V87KNw',
                      'm2-dC5QX8', 'm2-G68KXW', 'm2-QTwV7a', 'm2-pDAji8', 'm2-oxY6PQ', 'm2-AXsNSL', 'm2-nmHdMk',
                      'm2-mkGmpe', 'm2-gxDNRv', 'm2-oUT7WT', 'm2-MvvBk9', 'm2-i5wLSU', 'm2-s96RDD', 'm2-NdSZao',
                      'm2-gycK6z', 'm2-DhoSq8', 'm2-eyf48e', 'm2-9qnCXK', 'm2-V9rHyG', 'm2-5YkKkJ', 'm2-XP6vJf',
                      'm2-QEotEA', 'm2-7xxDyW', 'm2-37cQGL', 'm2-3qdHwU', 'm2-ZNzLxp', 'm2-YRBBfu', 'm2-yQUmMG',
                      'm2-y9sf9q', 'm2-uBCY5S', 'm2-TzR6v7', 'm2-QtXNo7', 'm2-QGowTB', 'm2-ozKj7w', 'm2-npCqtD',
                      'm2-mGefQv', 'm2-gyQtpZ', 'm2-G4bQxE', 'm2-xvyMxE', 'm2-Mj6hgu', 'm2-3KZQzJ', 'm2-cGQaRk',
                      'm2-bkrZsY', 'm2-BJn2VZ', 'm2-BF2ufp', 'm2-3SXzSV', 'm2-4kXick', 'm2-38nVTV', 'm2-3Ft7q9',
                      'm2-fChozQ', 'm2-n45rMk', 'm2-3UqUSN', 'm2-n4CyR5', 'm2-YztDkU', 'm2-9KEquH', 'm2-L3ukKo',
                      'm2-oDzvGR', 'm2-SccDgj', 'm2-2NXwCk', 'm2-oUUEge', 'm2-h54gk9', 'm2-TcfHz7', 'm2-i7rgAJ',
                      'm2-nHuJXs', 'm2-tDuWkr', 'm2-X4aXWG', 'm2-oD3bo4', 'm2-kLPe6f', 'm2-e3Nkuo', 'm2-ED6Fyo',
                      'm2-v72khW', 'm2-sSV2Bt', 'm2-pCrZWV', 'm2-zhhg9e', 'm2-53qCwY', 'm2-GRjNzo', 'm2-BvdErP',
                      'm2-34KSAo', 'm2-weJAnw', 'm2-y58Q75', 'm2-q8yvBm', 'm2-amRv2M', 'm2-3CuPYs', 'm2-ToGSRe',
                      'm2-jM3n9w', 'm2-Ho9wQu', 'm2-aEURHb', 'm2-ABxNA9', 'm2-ujTsiu', 'm2-NaDTTs', 'm2-MUZ8Cb',
                      'm2-jCKGBv', 'm2-39xr8K', 'm2-JBeRmh', 'm2-VnXdYZ', 'm2-mzGTDC', 'm2-xwByzt', 'm2-gEpiDa',
                      'm2-gvJsdn', 'm2-EYQGUT', 'm2-wgoGjP', 'm2-yEmM7C', 'm2-BP3DjK', 'm2-qFQ9kJ', 'm2-PQ5aTg',
                      'm2-ar4p3m', 'm2-Kouapg', 'm2-gxCe5X', 'm2-RgZWrD', 'm2-Uph9U3', 'm2-TUofv5', 'm2-ENT3qs',
                      'm2-MjqpEV', 'm2-Kyii3a', 'm2-LMD6ho', 'm2-HsBVJN', 'm2-Kr6366', 'm2-JqBYLh', 'm2-W4ZnHe',
                      'm2-yR5Hef', 'm2-ozGf4E', 'm2-CgiSAY', 'm2-BYjiah', 'm2-wF4k2T', 'm2-V8dKQC', 'm2-PVBn6f',
                      'm2-Mp25Fk', 'm2-MM8MkA', 'm2-2kH55m', 'm2-U2wv8C', 'm2-iZ7TLc', 'm2-2MHu4K', 'm2-m5XiUJ',
                      'm2-LX4eJE', 'm2-JiBCDX', 'm2-dLeBHd', 'm2-NqUrMx', 'm2-844sRF', 'm2-WydUVb', 'm2-YB7Zxj',
                      'm2-SrvHnQ', 'm2-pdN7uy', 'm2-6qaV6z', 'm2-ZNTJLD', 'm2-jwkH9W', 'm2-DBKZHX', 'm2-BbET5x',
                      'm2-7oH73S', 'm2-4Q8Stt', 'm2-Di2qm7', 'm2-QPbjvy', 'm2-boLTAg', 'm2-ND22rx', 'm2-FcN3wC',
                      'm2-RVmPKr', 'm2-f6biAx', 'm2-NESNdH', 'm2-8JqGbD', 'm2-sZ4qK3', 'm2-RCdbJH', 'm2-Y4wzR3',
                      'm2-nf9ieP', 'm2-AnXy9n', 'm2-MQhF9c', 'm2-yLHwrL', 'm2-c4imMZ', 'm2-G86emu', 'm2-5pS65p',
                      'm2-xw6CDX', 'm2-XMdsys', 'm2-tDscqy', 'm2-hHKDyB', 'm2-sVxWai', 'm2-Qb4JVo', 'm2-ahT4kF',
                      'm2-G8vc2N', 'm2-rhiSHJ', 'm2-BP5QSp', 'm2-Lq3fms', 'm2-voFSf6', 'm2-M5iqYC', 'm2-KtHTz8',
                      'm2-hjq9oe', 'm2-dgZP3L', 'm2-Q8Ba9A', 'm2-HAAUgG', 'm2-sJs93G', 'm2-GpvELC', 'm2-xxAxC9',
                      'm2-xMNGcb', 'm2-t9UZtV', 'm2-e3y7tj', 'm2-jQNpYD', 'm2-3URpmc', 'm2-4JSa55', 'm2-3fLBSP',
                      'm2-BdtP8F', 'm2-HXGr8V', 'm2-2YwKxe', 'm2-Xx7sqD', 'm2-3aogRq', 'm2-33eCP4', 'm2-hXgLde',
                      'm2-3wTS9S', 'm2-4Rt9sw', 'm2-3fGJE8', 'm2-8ERNNJ', 'm2-3TMBNi', 'm2-4SWHMc', 'm2-4Die3B',
                      'm2-3MgVni', 'm2-49m5VA', 'm2-E6qgW7', 'm2-X3sxhS', 'm2-4i2tZF', 'm2-zG7RFd', 'm2-MT4r8p',
                      'm2-4PF58Q', 'm2-36z874', 'm2-XisgNp', 'm2-aWAj59', 'm2-AWBv2b', 'm2-3Fa9rz', 'm2-hU5TjU',
                      'm2-Zu3u7N', 'm2-3MV28K', 'm2-3xWWe7', 'm2-9jS4Ro', 'm2-3H2Wig', 'm2-32supx', 'm2-YAdfpF',
                      'm2-4Lx8Lv', 'm2-3d2qwx', 'm2-gYFizT', 'm2-468ARL', 'm2-4HbGG2', 'm2-4R8qnW', 'm2-4W9wb5',
                      'm2-4VxKmQ', 'm2-3zPstM', 'v72khW', 'CcozBt', 'VduGeq', 's96RDD', 'kLPe6f', 'oUUEge', 'zhhg9e',
                      'sSV2Bt', 'NdSZao', 'wM9MFc', 'fChozQ', 'fngLcf', 'h54gk9', '3UqUSN', 'ujTsiu', 'icabqG',
                      'F34WtX', 'TcfHz7', 'eyf48e', 'HUBpy4', 'aEURHb', 'jCKGBv', 'BvdErP', 'jM3n9w', 'nHuJXs',
                      'pCrZWV', 'AnXy9n', '6RPWyT', 'qFQ9kJ', 'GbMzJy', 'yebiMQ', 'G86emu', 'sVxWai', 'hHKDyB',
                      'RCdbJH', 'Qf8cLC', 'SrvHnQ', 'pdN7uy', 'v3KDUF', 'X3sxhS', 'MT4r8p', '5fSMj3', 'pzbDPq',
                      '8ERNNJ', 'Kr6366', 'ycPovf', '8YYVHP', 'TCsq7r', 'T5VXtD', 'zxhtR4', 'qtZxkz', 'kBgpZ4',
                      'Sm9GYA', '9HhvMh', 'xrRPFj', 'VA9Yi5', '4rvG7f', 'xYnaUy', 'mzwDVC', 'SQZYHn', 'GpvELC',
                      'oUT7WT', 'g33df2', 'i5wLSU', 'opfXjQ', 'DwiwwU', 'ydSpNQ', '39kc9U', 'xHGR9u', 'oLXZ9a',
                      't9UZtV', 'q2oXnd', 'cyCNuA', '4SRYNi', '3SXzSV', '32WNFL', '3fLBSP', 'm6243p', '5cJWyf',
                      'WydUVb', '2MHu4K', 'TPgies', 'iZ7TLc', 'Ki8Aa3', 'Xo2ZNk', 'sFFUfJ', 'xwByzt', 'V8dKQC',
                      'V9rHyG', '6z2wf2', '2kH55m', 'LX4eJE', 'CdyN3k', 'yR5Hef', 'MM8MkA', 'dLeBHd', 'STGkhE',
                      'ozGf4E', 'ZNTJLD', 'Di2qm7', 'FcN3wC', 'uSTYZ2', '2hKeW2', 'FXFt6J', 'b99ooj', 'sAyJww',
                      'RVmPKr', 'gaQkyT', 'CgiSAY', '5YkKkJ', 'wF4k2T', 'yhE88h', 'DBKZHX', 'afGSAT', 'j3ZXpL',
                      'jwkH9W', 'ex4d6T', 'CLwScw', 'XisgNp', 'G4bQxE', 'jtjCjd', 'Zu3u7N', '4i2tZF', '4UF9K6',
                      'HXGr8V', 'q4bsMd', 'RyXauM', 'itJLqU', 'xwkVbK', 'tUrV5z', 'QGowTB', 'kwbwf8', 'npCqtD',
                      'xPAa2d', 'y9sf9q', 'E6qgW7', 'ozKj7w', 'yQUmMG', 'uBCY5S', 'HAAUgG', 'f8rfBe', 'NqUrMx',
                      'TzR6v7', 'Mp25Fk', '3H4uLZ', 'NMp4NB', 'x6NRYc', 'm5XiUJ', 'ND22rx', 'gyQtpZ', 'QtXNo7',
                      'BnbHLr', 'BJn2VZ', 'uTpHe5', 'KaJ5Pe', 'xvyMxE', 'dgZP3L', 'Ho9wQu', '8QXJTp', 'eAr6Lo',
                      'gkcF5C', 'MvvBk9', 'KtHTz8', 'k9vgiw', '4t7Ncr', 'iDqyV9', 'G8vc2N', 'mmYEfJ', 'hzNQRG',
                      'BP5QSp', 'ToGSRe', 'MUZ8Cb', 'gycK6z', '34KSAo', 'wVBQbT', 'cGQaRk', '3Ft7q9', 'zGrf4n',
                      '3eGjGj', '4JSa55', 'bkrZsY', '4SWHMc', '375EUk', '3Vw5ij', 'sZ4qK3', 'VT3bcF', 'DhoSq8',
                      '49m5VA', '5xXjjp', '3fGJE8', '33eCP4', 'YRBBfu', 'jQNpYD', 'V5iihT', 'sJs93G', 'ss4Fh8',
                      'SKJgwA', 'K2udH5', 'h2ZJd5', 'hjq9oe', 'Q8Ba9A', 'nf9ieP', '7xxDyW', 'sr4bvM', 'LosGHn',
                      'M5iqYC', 'ahT4kF', 'LCBdav', 'SpXNmN', 'Lq3fms', 'juVwzj', '3VC6ft', 'DxMevL', '3edUBd',
                      '4H6FKQ', '38nVTV']

        for sample_id in sample_ids:
            try:
                # find the sample by pseudonymized_id
                sample = Sample.objects.get(pseudonymized_id=sample_id)
            except Sample.DoesNotExist:
                print(f"{sample_id} ---> SAMPLE NOT FOUND")
                continue

            # find metadata linked to this sample
            metadata = Metadata.objects.filter(sample=sample).first()

            if metadata is None:
                print(f"{sample_id} ---> NO METADATA")
                continue

            # ent_date may be null
            if metadata.ent_date is None:
                print(f"{sample_id} ---> ent_date IS NULL")
            else:
                print(f"{sample_id} ---> {metadata.ent_date}")
